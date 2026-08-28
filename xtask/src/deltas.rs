//! The delta register, cross-checked between the source and the book.
//!
//! `book/src/deltas.md` is prose: what the original Fortran does, what this
//! port does instead, and the measured effect. The source carries the other
//! half, a `@delta` marker sitting on the code that actually deviates:
//!
//! ```text
//! // @port  teprob.f:1415-1442
//! // @delta D-001 class=B teprob.f:1439-1440
//! pub fn temperature_from_enthalpy(/* ... */) { /* ... */ }
//! ```
//!
//! Two halves that are maintained separately drift, and both failure directions
//! are silent. A delta documented with no marker reads as decided when nobody
//! can point at the line it applies to. A marker with no entry is a deviation
//! from the Fortran that never got its class, its measurement, or its sign-off.
//! So this module refuses both, by name, rather than reporting a count.
//!
//! It also emits `book/src/validation/delta-index.md`, which is the cheap half:
//! a table of every marker and where it sits. The table is convenient. The
//! cross-check is the point.

use std::collections::BTreeMap;
use std::fmt::Write as _;
use std::fs;
use std::path::Path;

use crate::report;
use crate::{LineRange, PROVENANCE_FILES, rust_sources, strip_comment_prefix, take_usize};

/// The marker that opens a delta claim, anchored at the start of a comment for
/// the same reason `@port` is: prose about the convention must not count.
const DELTA_MARKER: &str = "@delta";

/// The hand-written register the markers are checked against.
const REGISTER: &str = "book/src/deltas.md";

/// Where the generated index is written.
const INDEX: &str = "book/src/validation/delta-index.md";

/// A `@delta` marker found in the source.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Marker {
    /// The register id, normalised to `D-001` form.
    pub(crate) id: String,
    /// `A`, `B` or `C`, as `PLAN.org` defines them.
    pub(crate) class: char,
    /// Index into [`PROVENANCE_FILES`]: which vendored file the range is in.
    pub(crate) file: usize,
    pub(crate) range: LineRange,
}

/// A marker whose text is recognisably a claim but does not parse.
///
/// Reported rather than dropped: a mistyped marker is a marker that stopped
/// counting, and the failure mode is a delta that silently looks unmarked.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct MalformedMarker {
    pub(crate) text: String,
    pub(crate) why: &'static str,
}

/// One `## D-0NN` entry in `book/src/deltas.md`.
#[derive(Clone, Debug, PartialEq, Eq)]
pub(crate) struct Entry {
    pub(crate) id: String,
    /// The heading text with the id stripped off, for the index table.
    pub(crate) title: String,
    /// `A`, `B` or `C`, read from the `**Class X.**` line under the heading.
    pub(crate) class: char,
    /// One-based line of the `##` heading, so a failure can be navigated to.
    pub(crate) line: usize,
}

/// Where a marker was found, for the index and for failure messages.
#[derive(Clone, Debug)]
pub(crate) struct Located {
    pub(crate) marker: Marker,
    /// Workspace-relative path of the Rust file.
    pub(crate) path: String,
    pub(crate) line: usize,
}

pub(crate) fn cmd_deltas(root: &Path) -> Result<(), String> {
    let found = scan_markers(root)?;
    let register = root.join(REGISTER);
    let text = fs::read_to_string(&register).map_err(|e| {
        format!(
            "reading {REGISTER}: {e}\nThe delta index is generated against the \
             register, so without it there is nothing to check."
        )
    })?;
    let entries = parse_register(&text)?;

    // Write the page before checking, so a failing cross-check still leaves a
    // readable artefact naming what it found. The check is the deliverable, but
    // a tool that reports a problem and then destroys the evidence is annoying.
    let page = render_index(root, &entries, &found);
    report::write_generated(root, INDEX, &page)?;

    cross_check(&entries, &found)?;

    println!(
        "[ok] deltas: {} register entr(ies), {} marker(s), classes agree",
        entries.len(),
        found.len()
    );
    Ok(())
}

/// Every `@delta` marker in the porting crates.
pub(crate) fn scan_markers(root: &Path) -> Result<Vec<Located>, String> {
    let mut out = Vec::new();
    let mut malformed = Vec::new();
    for file in rust_sources(root)? {
        let text =
            fs::read_to_string(&file).map_err(|e| format!("reading {}: {e}", file.display()))?;
        let relative = relative(root, &file);
        for (index, line) in text.lines().enumerate() {
            match parse_marker_line(line) {
                Some(Ok(marker)) => out.push(Located {
                    marker,
                    path: relative.clone(),
                    line: index + 1,
                }),
                Some(Err(bad)) => malformed.push((relative.clone(), index + 1, bad)),
                None => {}
            }
        }
    }

    if !malformed.is_empty() {
        let lines: Vec<String> = malformed
            .iter()
            .map(|(path, line, bad)| format!("  {path}:{line}: {} ({})", bad.text.trim(), bad.why))
            .collect();
        return Err(format!(
            "{} malformed `@delta` marker(s):\n{}\n\
             The shape is `// @delta D-001 class=B teprob.f:1439-1440`. A marker \
             that does not parse is a marker that stopped counting, and the \
             delta it belongs to then reads as undocumented in the source.",
            malformed.len(),
            lines.join("\n")
        ));
    }

    out.sort_by(|a, b| {
        a.marker
            .id
            .cmp(&b.marker.id)
            .then_with(|| a.path.cmp(&b.path))
            .then_with(|| a.line.cmp(&b.line))
    });
    Ok(out)
}

/// Parse one line as `// @delta D-001 class=B teprob.f:1439-1440`.
///
/// `None` for anything that is not a marker at all. `Some(Err(..))` for a line
/// that opens with the marker and then does not parse, which is the case worth
/// shouting about.
pub(crate) fn parse_marker_line(line: &str) -> Option<Result<Marker, MalformedMarker>> {
    let body = strip_comment_prefix(line.trim_start())?;
    let after = body.trim_start().strip_prefix(DELTA_MARKER)?;
    // Whitespace after the marker, so `@deltas` is not a marker.
    if !after.starts_with(char::is_whitespace) {
        return None;
    }
    let text = after.trim().to_string();
    let bad = |why| {
        Some(Err(MalformedMarker {
            text: text.clone(),
            why,
        }))
    };

    let mut words = text.split_whitespace();
    let Some(id) = words.next().and_then(normalise_id) else {
        return bad("the id must look like `D-001`");
    };
    let Some(class) = words.next().and_then(parse_class) else {
        return bad("expected `class=A`, `class=B` or `class=C`");
    };
    let Some(rest) = words.next() else {
        return bad("expected a vendored source range, e.g. `teprob.f:1439-1440`");
    };
    let Some((file, tail)) = PROVENANCE_FILES
        .iter()
        .enumerate()
        .find_map(|(index, (tag, _))| rest.strip_prefix(tag).map(|r| (index, r)))
    else {
        return bad("the source range must name a vendored Fortran file");
    };
    let Some((start, tail)) = take_usize(tail) else {
        return bad("the source range must start with a line number");
    };
    let end = tail
        .strip_prefix('-')
        .and_then(take_usize)
        .map_or(start, |(e, _)| e);

    Some(Ok(Marker {
        id,
        class,
        file,
        range: LineRange {
            start,
            end: end.max(start),
        },
    }))
}

/// `D-1`, `D-01` and `D-001` are the same delta; the register writes `D-001`.
///
/// Normalising rather than demanding one spelling means a marker written
/// `D-11` still matches its entry instead of being reported as an orphan,
/// which would be a confusing failure for a cosmetic mistake.
fn normalise_id(word: &str) -> Option<String> {
    let digits = word.strip_prefix("D-")?;
    if digits.is_empty() || !digits.bytes().all(|b| b.is_ascii_digit()) {
        return None;
    }
    let number: u32 = digits.parse().ok()?;
    Some(format!("D-{number:03}"))
}

fn parse_class(word: &str) -> Option<char> {
    let value = word.strip_prefix("class=")?;
    let mut chars = value.chars();
    let class = chars.next()?;
    // Exactly one letter: `class=BB` is a typo, not a class.
    (chars.next().is_none() && matches!(class, 'A' | 'B' | 'C')).then_some(class)
}

/// The `## D-0NN` entries of `book/src/deltas.md`, with their classes.
///
/// The class is read from the `**Class X.**` line under the heading rather than
/// from the heading itself, because that line is where the register actually
/// states it and a second copy in the heading would be one more thing to drift.
pub(crate) fn parse_register(text: &str) -> Result<Vec<Entry>, String> {
    let mut out: Vec<Entry> = Vec::new();
    let mut pending: Option<(String, String, usize)> = None;
    let mut fenced = false;

    for (index, raw) in text.lines().enumerate() {
        let line = raw.trim();
        // The register documents the marker convention with a fenced example.
        // Headings and class lines inside a fence are illustration, not entries.
        if line.starts_with("```") {
            fenced = !fenced;
            continue;
        }
        if fenced {
            continue;
        }

        if let Some(rest) = line.strip_prefix("## ") {
            if let Some((id, title)) = split_heading(rest) {
                if let Some((id, _, line_no)) = pending.take() {
                    return Err(missing_class(&id, line_no));
                }
                pending = Some((id, title, index + 1));
            }
            continue;
        }

        if let Some((id, title, line_no)) = pending.take() {
            match parse_class_line(line) {
                Some(class) => out.push(Entry {
                    id,
                    title,
                    class,
                    line: line_no,
                }),
                // Blank lines between the heading and the class line are normal.
                None if line.is_empty() => pending = Some((id, title, line_no)),
                None => return Err(missing_class(&id, line_no)),
            }
        }
    }
    if let Some((id, _, line_no)) = pending {
        return Err(missing_class(&id, line_no));
    }

    // Two entries with one id would make the cross-check ambiguous, and the
    // second one's class would silently win.
    let mut seen: BTreeMap<&str, usize> = BTreeMap::new();
    for entry in &out {
        if let Some(first) = seen.insert(&entry.id, entry.line) {
            return Err(format!(
                "{REGISTER} has two entries for {} (lines {first} and {}). \
                 One id, one entry.",
                entry.id, entry.line
            ));
        }
    }
    Ok(out)
}

fn missing_class(id: &str, line: usize) -> String {
    format!(
        "{REGISTER}:{line}: the entry for {id} has no `**Class X.**` line under \
         its heading.\nEvery entry states its class there, because the class is \
         what decides whether the deviation could ship as the default."
    )
}

/// `D-001 — <title>` from a `## ` heading, or `None` for any other heading.
fn split_heading(rest: &str) -> Option<(String, String)> {
    let mut words = rest.split_whitespace();
    let id = normalise_id(words.next()?)?;
    let title = rest
        .split_once(char::is_whitespace)
        .map_or("", |(_, tail)| tail)
        .trim_start_matches(['\u{2014}', '\u{2013}', '-', ' '])
        .trim()
        .to_string();
    Some((id, title))
}

/// The class from a `**Class B.** ...` line.
///
/// Tolerant on purpose: entries write `**Class A, and load-bearing.**` and
/// `**Class C. Reproduced by default; ...`, and the letter is the part that
/// matters. Requiring an exact shape would turn a rewording into a gate failure.
fn parse_class_line(line: &str) -> Option<char> {
    let rest = line.strip_prefix("**Class ")?;
    let mut chars = rest.chars();
    let class = chars.next()?;
    if !matches!(class, 'A' | 'B' | 'C') {
        return None;
    }
    // A word boundary, so `**Class Bx**` is not read as class B.
    match chars.next() {
        Some('.') | Some(',') | Some(' ') | Some('*') | None => Some(class),
        _ => None,
    }
}

/// The whole point: neither half of the register may drift from the other.
fn cross_check(entries: &[Entry], found: &[Located]) -> Result<(), String> {
    let documented: BTreeMap<&str, &Entry> = entries.iter().map(|e| (e.id.as_str(), e)).collect();
    let mut marked: BTreeMap<&str, Vec<&Located>> = BTreeMap::new();
    for located in found {
        marked
            .entry(located.marker.id.as_str())
            .or_default()
            .push(located);
    }

    let mut problems = Vec::new();

    for (id, entry) in &documented {
        if !marked.contains_key(id) {
            problems.push(format!(
                "  {id} is documented at {REGISTER}:{} but no `@delta` marker \
                 names it.\n    Nothing in the source points at the code that \
                 deviates, so the entry cannot be\n    checked against what the \
                 port does, and deleting the deviation would not\n    disturb \
                 the register at all. Put `// @delta {id} class={} \
                 <file>:<lines>`\n    on the item the entry is about.",
                entry.line, entry.class
            ));
        }
    }

    for (id, sites) in &marked {
        let Some(entry) = documented.get(id) else {
            let where_at: Vec<String> = sites
                .iter()
                .map(|s| format!("{}:{}", s.path, s.line))
                .collect();
            problems.push(format!(
                "  {id} is marked at {} but {REGISTER} has no `## {id}` \
                 entry.\n    A deviation from the Fortran with no register entry \
                 has no class, no measured\n    effect and no sign-off, which for \
                 a Class C deviation is the difference between\n    reproducing \
                 the benchmark and quietly changing it.",
                where_at.join(", ")
            ));
            continue;
        };
        for site in sites {
            if site.marker.class != entry.class {
                problems.push(format!(
                    "  {id} is class {} at {}:{} but class {} at {REGISTER}:{}.\n\
                     \x20   The class decides the disposition, so the two \
                     disagreeing means one of them is\n    granting a deviation \
                     a licence the other refuses it.",
                    site.marker.class, site.path, site.line, entry.class, entry.line
                ));
            }
        }
    }

    if problems.is_empty() {
        return Ok(());
    }
    Err(format!(
        "the delta register and the source disagree in {} place(s):\n\n{}\n\n\
         Both halves are needed: the register says what the deviation is and \
         what it costs,\nthe marker says which line it is. Either alone rots.",
        problems.len(),
        problems.join("\n\n")
    ))
}

fn render_index(root: &Path, entries: &[Entry], found: &[Located]) -> String {
    let mut marked: BTreeMap<&str, Vec<&Located>> = BTreeMap::new();
    for located in found {
        marked
            .entry(located.marker.id.as_str())
            .or_default()
            .push(located);
    }
    let mut ids: Vec<&str> = entries.iter().map(|e| e.id.as_str()).collect();
    for id in marked.keys() {
        if !ids.contains(id) {
            ids.push(id);
        }
    }
    ids.sort_unstable();

    let mut page = report::header(
        root,
        "Delta marker index",
        "cargo xtask deltas",
        "Every row was collected from the source and from the register by that \
         run.",
    );
    page.push_str(
        "\nEvery deliberate deviation from the original Fortran carries two \
         things: an entry\nin the [quirk and delta register](../deltas.md), and \
         a `@delta` marker on the code\nit applies to. This table is collected \
         from the markers and checked against the\nregister. `cargo xtask \
         deltas` fails if an entry has no marker, if a marker has no\nentry, or \
         if the two disagree about the class.\n\n",
    );

    let by_class = |wanted: char| entries.iter().filter(|e| e.class == wanted).count();
    let _ = writeln!(
        page,
        "{} register entries, {} markers across {} source location(s): \
         {} class A, {} class B, {} class C.\n",
        entries.len(),
        found.len(),
        marked.values().map(Vec::len).sum::<usize>(),
        by_class('A'),
        by_class('B'),
        by_class('C'),
    );

    page.push_str("| delta | class | deviates from | marked at | register |\n");
    page.push_str("|---|---|---|---|---|\n");
    for id in ids {
        let entry = entries.iter().find(|e| e.id == id);
        let sites = marked.get(id).cloned().unwrap_or_default();
        let class = entry
            .map(|e| e.class.to_string())
            .unwrap_or_else(|| "**undocumented**".to_string());
        let fortran = if sites.is_empty() {
            "**unmarked**".to_string()
        } else {
            let mut ranges: Vec<String> = sites
                .iter()
                .map(|s| {
                    format!(
                        "`{}:{}`",
                        PROVENANCE_FILES[s.marker.file].0.trim_end_matches(':'),
                        s.marker.range
                    )
                })
                .collect();
            ranges.dedup();
            ranges.join("<br>")
        };
        let at = if sites.is_empty() {
            "-".to_string()
        } else {
            sites
                .iter()
                .map(|s| format!("`{}:{}`", s.path, s.line))
                .collect::<Vec<_>>()
                .join("<br>")
        };
        let title = entry.map_or("*no register entry*", |e| e.title.as_str());
        let _ = writeln!(page, "| {id} | {class} | {fortran} | {at} | {title} |");
    }

    page.push_str(
        "\nThe class column comes from the register. A row reading \
         **unmarked** or\n**undocumented** is not a formatting artefact: it is \
         a half of the register that is\nmissing, and the run that wrote this \
         page exited non-zero saying so. The page is\nwritten before the \
         cross-check so that the evidence survives the failure.\n\n\
         Class A has no numerical effect, class B is numerically observable and \
         fixed with a\nmeasured delta, and class C is behaviour-defining, \
         reproduced by default and changed\nonly behind a flag on explicit \
         sign-off. `PLAN.org` defines them.\n",
    );
    page
}

fn relative(root: &Path, path: &Path) -> String {
    path.strip_prefix(root)
        .unwrap_or(path)
        .to_string_lossy()
        .replace('\\', "/")
}

#[cfg(test)]
mod tests {
    use super::*;

    fn marker(line: &str) -> Marker {
        parse_marker_line(line)
            .expect("recognised as a marker")
            .expect("well formed")
    }

    #[test]
    fn parses_a_marker() {
        let m = marker("// @delta D-001 class=B teprob.f:1439-1440");
        assert_eq!(m.id, "D-001");
        assert_eq!(m.class, 'B');
        assert_eq!(m.file, 0);
        assert_eq!(
            m.range,
            LineRange {
                start: 1439,
                end: 1440
            }
        );
    }

    #[test]
    fn parses_an_indented_marker_naming_the_driver() {
        let m = marker("        // @delta D-011 class=C temain_mod.f:366-368");
        assert_eq!(m.id, "D-011");
        assert_eq!(m.class, 'C');
        assert_eq!(m.file, 1);
    }

    #[test]
    fn a_single_line_range_is_a_range_of_one() {
        let m = marker("// @delta D-006 class=A teprob.f:698");
        assert_eq!(
            m.range,
            LineRange {
                start: 698,
                end: 698
            }
        );
    }

    #[test]
    fn ids_normalise_so_a_cosmetic_spelling_still_matches() {
        assert_eq!(marker("// @delta D-1 class=A teprob.f:1").id, "D-001");
        assert_eq!(marker("// @delta D-01 class=A teprob.f:1").id, "D-001");
    }

    /// Prose about the convention must not be collected, the same false
    /// positive `@port` had.
    #[test]
    fn prose_and_glued_words_are_not_markers() {
        assert!(parse_marker_line("//! markers look like `@delta D-001 ...`").is_none());
        assert!(parse_marker_line("// @deltas D-001 class=A teprob.f:1").is_none());
        assert!(parse_marker_line("let x = 1; // @delta").is_none());
        assert!(parse_marker_line("// see D-001").is_none());
    }

    /// A marker that opens correctly and then does not parse is loud, not
    /// silent: silently dropping it makes the delta look unmarked.
    #[test]
    fn a_malformed_marker_is_reported_rather_than_dropped() {
        for (line, fragment) in [
            ("// @delta D001 class=A teprob.f:1", "id"),
            ("// @delta D-001 clas=A teprob.f:1", "class="),
            ("// @delta D-001 class=Z teprob.f:1", "class="),
            ("// @delta D-001 class=A", "vendored source range"),
            (
                "// @delta D-001 class=A teprob.for:1",
                "vendored Fortran file",
            ),
            ("// @delta D-001 class=A teprob.f:x", "line number"),
        ] {
            let parsed = parse_marker_line(line).expect("recognised as a marker");
            let error = parsed.expect_err(line);
            assert!(error.why.contains(fragment), "{line}: {}", error.why);
        }
    }

    #[test]
    fn reads_entries_and_their_classes() {
        let text = "# Register\n\n\
                    ## D-001 \u{2014} `TESUB2` reports success\n\n\
                    **Class B.** `teprob.f:1439-1440`.\n\n\
                    ### What the original does\n\n\
                    ## D-010 \u{2014} stale measurements\n\n\
                    **Class A, and load-bearing.** `temain_mod.f:366-411`.\n";
        let entries = parse_register(text).expect("parses");
        assert_eq!(entries.len(), 2);
        assert_eq!(entries[0].id, "D-001");
        assert_eq!(entries[0].class, 'B');
        assert_eq!(entries[0].title, "`TESUB2` reports success");
        assert_eq!(entries[1].id, "D-010");
        assert_eq!(entries[1].class, 'A');
    }

    /// The register shows the convention in a fenced example. Counting that as
    /// an entry would invent a D-001 with no heading and no class.
    #[test]
    fn a_fenced_example_is_not_an_entry() {
        let text = "```rust\n\
                    // @delta D-001 class=B teprob.f:1439-1440\n\
                    ## D-999 \u{2014} not an entry\n\
                    **Class A.** nope\n\
                    ```\n\
                    ## D-002 \u{2014} real\n\n**Class A.** `teprob.f:1`.\n";
        let entries = parse_register(text).expect("parses");
        assert_eq!(entries.len(), 1);
        assert_eq!(entries[0].id, "D-002");
    }

    #[test]
    fn an_entry_with_no_class_line_is_an_error() {
        let text = "## D-001 \u{2014} a title\n\n### What the original does\n";
        let error = parse_register(text).expect_err("no class");
        assert!(error.contains("no `**Class X.**`"), "{error}");
    }

    #[test]
    fn two_entries_with_one_id_are_an_error() {
        let text =
            "## D-001 \u{2014} a\n\n**Class A.** x\n\n## D-001 \u{2014} b\n\n**Class A.** y\n";
        let error = parse_register(text).expect_err("duplicate");
        assert!(error.contains("two entries"), "{error}");
    }

    fn located(id: &str, class: char) -> Located {
        Located {
            marker: Marker {
                id: id.to_string(),
                class,
                file: 0,
                range: LineRange { start: 1, end: 1 },
            },
            path: "crates/tepsim-core/src/thermo.rs".to_string(),
            line: 370,
        }
    }

    fn entry(id: &str, class: char) -> Entry {
        Entry {
            id: id.to_string(),
            title: "a deviation".to_string(),
            class,
            line: 46,
        }
    }

    #[test]
    fn matching_halves_pass() {
        cross_check(&[entry("D-001", 'B')], &[located("D-001", 'B')]).expect("agree");
    }

    /// The two failures this command exists for.
    #[test]
    fn documented_but_unmarked_fails_by_name() {
        let error = cross_check(&[entry("D-008", 'A')], &[]).expect_err("no marker");
        assert!(error.contains("D-008"), "{error}");
        assert!(error.contains("no `@delta` marker"), "{error}");
    }

    #[test]
    fn marked_but_undocumented_fails_by_name() {
        let error = cross_check(&[], &[located("D-042", 'C')]).expect_err("no entry");
        assert!(error.contains("D-042"), "{error}");
        assert!(error.contains("no `## D-042` entry"), "{error}");
        // And it names the file, because the fix is to write the entry and the
        // author needs to see what the marker is sitting on.
        assert!(error.contains("thermo.rs:370"), "{error}");
    }

    #[test]
    fn a_class_disagreement_fails() {
        let error =
            cross_check(&[entry("D-001", 'A')], &[located("D-001", 'B')]).expect_err("classes");
        assert!(error.contains("class B at"), "{error}");
        assert!(error.contains("class A at"), "{error}");
    }

    #[test]
    fn several_problems_are_all_reported_not_just_the_first() {
        let error = cross_check(
            &[entry("D-001", 'A'), entry("D-002", 'A')],
            &[located("D-003", 'C')],
        )
        .expect_err("three problems");
        for id in ["D-001", "D-002", "D-003"] {
            assert!(error.contains(id), "{id} missing from:\n{error}");
        }
        assert!(error.contains("3 place(s)"), "{error}");
    }

    #[test]
    fn class_lines_are_read_tolerantly_but_not_loosely() {
        assert_eq!(parse_class_line("**Class B.** `teprob.f:1`."), Some('B'));
        assert_eq!(
            parse_class_line("**Class A, and load-bearing.**"),
            Some('A')
        );
        assert_eq!(
            parse_class_line("**Class C. Reproduced by default"),
            Some('C')
        );
        assert_eq!(parse_class_line("**Class D.**"), None);
        assert_eq!(parse_class_line("**Classes A and B**"), None);
        assert_eq!(parse_class_line("The class is B"), None);
    }
}
