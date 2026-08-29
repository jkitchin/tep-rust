//! The book's code listings are the example files, and its notebook links go
//! to notebooks that exist and have been run.
//!
//! B-0068. `book/src/tutorials/scheduling-a-fault.md` used to carry a
//! `rust,ignore` block and a transcript pasted beside it, and both went stale
//! the moment a default changed: the transcript showed `active [4, 12]` from
//! the driver's forced `IDV(12)`, which stopped happening. Nothing could have
//! caught that, because nothing compiled the listing or ran it.
//!
//! B-0077a moved the worked examples out of the book altogether. The four
//! tutorial pages are narrative now and link to `notebooks/*.ipynb`, which are
//! executed, committed with their outputs and plots, and rendered beside the
//! book by `.github/workflows/pages.yml`. So this file lost three listings and
//! gained the two checks that failure mode needs instead: a link from a page to
//! a notebook that is not there, and a notebook committed without having been
//! run.
//!
//! One listing is left, the quickstart on `book/src/python.md`, and it is
//! pinned exactly as before. It is the only code in the book, and leaving it
//! unpinned in the same iteration that removed three listings *because* they
//! had gone stale would be a strange conclusion to draw.
//!
//! Two findings from the move, both of which this file would now catch.
//! `book/src/tutorials/first-run.md` was never pinned and quoted a scenario
//! digest of `b3415d9a395b8c70` for a run that had produced `4400cc4f5f4f570e`
//! ever since the `driver_forces_idv12` default flipped. And
//! `scheduling-a-fault.md`'s prose quoted a digest of `246b561b64789ada` that
//! did not match the transcript printed six lines above it.

use std::fs;
use std::path::{Path, PathBuf};
use std::process::Command;

fn repo(rest: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rest)
}

/// Where a page's listing really lives, and how to run it.
enum Listing {
    /// A compiled example, `crates/<krate>/examples/<name>.rs`. Shown in a
    /// ```` ```rust ```` fence and run by `cargo run --example`.
    Rust {
        krate: &'static str,
        name: &'static str,
    },
    /// A script under `book/examples/`. Shown in a ```` ```python ```` fence
    /// and run by the virtualenv `cargo xtask python` builds.
    Python { script: &'static str },
}

/// A page, the listing it shows, and whether it also quotes that listing's
/// output.
///
/// `README.md` shows a snippet and no transcript, so it is pinned for source
/// and skipped for output.
struct Listed {
    page: &'static str,
    listing: Listing,
    transcript: bool,
}

const LISTED: &[Listed] = &[
    Listed {
        page: "book/src/python.md",
        listing: Listing::Python {
            script: "book/examples/quickstart.py",
        },
        transcript: true,
    },
    Listed {
        page: "README.md",
        listing: Listing::Rust {
            krate: "tepsim",
            name: "readme_snippet",
        },
        transcript: false,
    },
];

/// The body of the first fence opening with `tag` on a page.
fn first_block(markdown: &str, tag: &str) -> String {
    let mut lines = markdown.lines();
    let mut body = Vec::new();
    for line in lines.by_ref() {
        if line.starts_with(tag) {
            break;
        }
    }
    for line in lines {
        if line.starts_with("```") {
            return body.join("\n");
        }
        body.push(line);
    }
    panic!("the fence was never closed");
}

/// One example crate's source path.
fn example_path(krate: &str, name: &str) -> PathBuf {
    repo(&format!("crates/{krate}/examples/{name}.rs"))
}

/// The body of `fn main`, dedented by one level.
///
/// A page that is teaching a whole program shows the `fn main` wrapper; a
/// quick-start snippet should not have to, and making `README.md` carry one
/// just to satisfy a test would be the test dictating the prose. So a page
/// whose fence has no `fn main` is compared against the body of the example's,
/// which is the same code either way.
fn main_body(source: &str) -> String {
    let Some(start) = source.find("fn main() {") else {
        return source.to_string();
    };
    let after = &source[start + "fn main() {".len()..];
    let Some(end) = after.rfind('}') else {
        return source.to_string();
    };
    let mut out: Vec<String> = Vec::new();
    for line in after[..end].lines() {
        out.push(line.strip_prefix("    ").unwrap_or(line).to_string());
    }
    // The `use` lines live above `fn main` and belong to the snippet too.
    let mut head: Vec<&str> = source[..start]
        .lines()
        .filter(|l| !l.starts_with("//!"))
        .collect();
    while head.first().is_some_and(|l| l.trim().is_empty()) {
        head.remove(0);
    }
    while out.first().is_some_and(|l| l.trim().is_empty()) {
        out.remove(0);
    }
    while out.last().is_some_and(|l| l.trim().is_empty()) {
        out.pop();
    }
    let head = head.join("\n");
    let head = head.trim_end();
    if head.is_empty() {
        out.join("\n")
    } else {
        // One blank line between the imports and the body, which is how the
        // source has it and how a snippet reads.
        format!("{head}\n\n{}", out.join("\n"))
    }
}

/// A Python script without the module docstring at the top of it.
///
/// The docstring says which page the file is pinned to and how to run it. That
/// is about the pinning rather than about the page and would be noise in the
/// book, so it is stripped here for exactly the reason the Rust examples' `//!`
/// headers are.
fn without_module_docstring(source: &str) -> String {
    let mut lines = source.lines();
    let Some(first) = lines.next() else {
        return String::new();
    };
    if !first.starts_with("\"\"\"") {
        return source.to_string();
    }
    // A one-line docstring opens and closes on the same line.
    if !(first.len() > 3 && first.ends_with("\"\"\"")) {
        for line in lines.by_ref() {
            if line.trim_end().ends_with("\"\"\"") {
                break;
            }
        }
    }
    lines
        .skip_while(|l| l.trim().is_empty())
        .collect::<Vec<_>>()
        .join("\n")
}

/// The listing as it is shown on the page, and as the file on disk has it.
fn shown_and_expected(listed: &Listed, markdown: &str) -> (String, String) {
    match &listed.listing {
        Listing::Rust { krate, name } => {
            let example = fs::read_to_string(example_path(krate, name))
                .unwrap_or_else(|_| panic!("crates/{krate}/examples/{name}.rs"));
            // The page omits the file's own `//!` header, which is about the
            // pinning rather than about the page.
            let body: String = example
                .lines()
                .skip_while(|l| l.starts_with("//!") || l.is_empty())
                .collect::<Vec<_>>()
                .join("\n");
            let shown = first_block(markdown, "```rust");
            // A page that shows no `fn main` is showing its body.
            let body = if shown.contains("fn main") {
                body
            } else {
                main_body(&body)
            };
            (shown, body)
        }
        Listing::Python { script } => {
            let source = fs::read_to_string(repo(script)).unwrap_or_else(|_| panic!("{script}"));
            (
                first_block(markdown, "```python"),
                without_module_docstring(&source),
            )
        }
    }
}

#[test]
fn every_listing_is_the_example_file() {
    for listed in LISTED {
        let page = listed.page;
        let markdown = fs::read_to_string(repo(page)).unwrap_or_else(|_| panic!("{page}"));
        let (shown, expected) = shown_and_expected(listed, &markdown);
        assert!(
            shown.lines().count() >= 4,
            "{page}: the parse found {} lines, so it is reading the wrong fence",
            shown.lines().count()
        );
        assert_eq!(
            shown.trim_end(),
            expected.trim_end(),
            "{page} and its example file have drifted"
        );
    }
}

/// Whether to run the expensive half: the listings, and the notebooks.
///
/// Off by default and on in `cargo xtask validate`, in the same spirit as
/// `xtask fidelity` versus the live oracle diff. The cheap checks that catch
/// gross breakage run on every commit and are unconditional, because a listing
/// drifting from its file and a link pointing at nothing are the failures that
/// actually happen.
fn expensive_checks_requested() -> bool {
    std::env::var("TEP_BOOK").is_ok_and(|v| v != "0")
}

/// Where `cargo xtask python` leaves the interpreter with `tepsim` in it.
#[cfg(windows)]
const VENV_PYTHON: &str = ".xtask-python/venv/Scripts/python.exe";
#[cfg(not(windows))]
const VENV_PYTHON: &str = ".xtask-python/venv/bin/python";

/// That interpreter, if it exists.
///
/// It usually does not, when `cargo test --workspace` runs: `xtask ci` builds
/// the wheel *after* the tests, so the venv the gate creates does not exist
/// while the gate's tests are running. A missing interpreter is therefore a
/// skip rather than a failure, on the same reasoning the oracle job already
/// uses for a missing gfortran. Run `cargo xtask python` and then
/// `cargo xtask validate`, in that order, to get the check.
fn venv_python() -> Option<PathBuf> {
    let path = repo(VENV_PYTHON);
    path.is_file().then_some(path)
}

/// And the transcript the page quotes is the output that listing produces.
///
/// Run rather than trusted. Without this the code could be pinned and the
/// numbers beside it still be from a build two defaults ago, which is exactly
/// what happened: `scheduling-a-fault.md` showed `active [4, 12]` from the
/// driver's forced `IDV(12)` for a day after that stopped being the default.
#[test]
fn every_transcript_is_the_output_it_quotes() {
    if !expensive_checks_requested() {
        println!("skipped: set TEP_BOOK=1 to re-run the listings; `cargo xtask validate` does");
        return;
    }
    let interpreter = venv_python();
    for listed in LISTED {
        if !listed.transcript {
            continue;
        }
        let page = listed.page;
        let markdown = fs::read_to_string(repo(page)).unwrap_or_else(|_| panic!("{page}"));
        let quoted = markdown
            .split_once("```text\n")
            .expect("a text fence")
            .1
            .split_once("```")
            .expect("a closing fence")
            .0;

        let output = match &listed.listing {
            Listing::Rust { krate, name } => Command::new(env!("CARGO"))
                .args([
                    "run",
                    "--quiet",
                    "--offline",
                    "-p",
                    krate,
                    "--example",
                    name,
                ])
                .current_dir(repo("."))
                .output()
                .expect("the example runs"),
            Listing::Python { script } => {
                let Some(python) = interpreter.as_ref() else {
                    println!(
                        "skipped {page}: {VENV_PYTHON} does not exist. \
                         `cargo xtask python` builds it; run that first."
                    );
                    continue;
                };
                Command::new(python)
                    .arg(repo(script))
                    .current_dir(repo("."))
                    .output()
                    .expect("the script runs")
            }
        };
        assert!(
            output.status.success(),
            "{page}'s example failed: {}",
            String::from_utf8_lossy(&output.stderr)
        );
        let produced = String::from_utf8(output.stdout).expect("utf-8");

        assert!(
            quoted.lines().count() > 15,
            "{page}: the quoted transcript is {} lines, so the parse is wrong",
            quoted.lines().count()
        );
        assert_eq!(
            produced.trim_end(),
            quoted.trim_end(),
            "{page} quotes output its example no longer produces"
        );
    }
}

// ---------------------------------------------------------------------------
// the notebooks the book sends readers to
// ---------------------------------------------------------------------------

/// Every markdown file under `book/src`, as its path relative to that root.
fn book_pages() -> Vec<PathBuf> {
    fn walk(dir: &Path, root: &Path, out: &mut Vec<PathBuf>) {
        let entries = fs::read_dir(dir).unwrap_or_else(|e| panic!("{}: {e}", dir.display()));
        for entry in entries {
            let path = entry.expect("a directory entry").path();
            if path.is_dir() {
                walk(&path, root, out);
            } else if path.extension().is_some_and(|e| e == "md") {
                out.push(
                    path.strip_prefix(root)
                        .expect("under the root")
                        .to_path_buf(),
                );
            }
        }
    }
    let root = repo("book/src");
    let mut out = Vec::new();
    walk(&root, &root, &mut out);
    out.sort();
    out
}

/// Every `](target)` on a page.
fn link_targets(markdown: &str) -> Vec<&str> {
    let mut out = Vec::new();
    let mut rest = markdown;
    while let Some(open) = rest.find("](") {
        let after = &rest[open + 2..];
        let Some(close) = after.find(')') else { break };
        out.push(&after[..close]);
        rest = &after[close..];
    }
    out
}

/// A link resolved against the page carrying it, as a site-relative path.
///
/// This is not pedantry. `book/src/python.md` renders to `<site>/python.html`
/// and `book/src/tutorials/a-detector.md` to `<site>/tutorials/a-detector.html`,
/// so the *same* notebook is `notebooks/x.html` from one page and
/// `../notebooks/x.html` from the other. Getting that wrong produces a link
/// that is present, plausible and 404, which no existence check on the notebook
/// alone would catch.
fn resolve(page: &Path, target: &str) -> Option<Vec<String>> {
    let mut path: Vec<String> = page
        .parent()?
        .components()
        .map(|c| c.as_os_str().to_string_lossy().into_owned())
        .collect();
    for part in target.split('/') {
        match part {
            "" | "." => {}
            ".." => {
                path.pop()?;
            }
            other => path.push(other.to_string()),
        }
    }
    Some(path)
}

/// Every notebook the book links to exists, at the depth the link claims.
#[test]
fn every_notebook_link_resolves() {
    let mut checked = 0;
    for page in book_pages() {
        let markdown = fs::read_to_string(repo("book/src").join(&page))
            .unwrap_or_else(|e| panic!("{}: {e}", page.display()));
        for target in link_targets(&markdown) {
            if !target.contains("notebooks/") || !target.ends_with(".html") {
                continue;
            }
            let resolved = resolve(&page, target).unwrap_or_else(|| {
                panic!(
                    "{} links to {target}, which escapes the book root",
                    page.display()
                )
            });
            assert_eq!(
                resolved.len(),
                2,
                "{} links to {target}, which resolves to {resolved:?}; the \
                 rendered notebooks live at <site>/notebooks/",
                page.display()
            );
            assert_eq!(
                resolved[0],
                "notebooks",
                "{} links to {target}, which resolves to {resolved:?}",
                page.display()
            );
            let notebook = resolved[1].replace(".html", ".ipynb");
            let path = repo("notebooks").join(&notebook);
            assert!(
                path.is_file(),
                "{} links to {target}, but notebooks/{notebook} does not exist",
                page.display()
            );
            checked += 1;
        }
    }
    assert!(
        checked >= 4,
        "found only {checked} notebook links in the book, so either the book \
         stopped linking to them or `link_targets` stopped finding them"
    );
}

/// Every notebook was executed before it was committed.
///
/// The notebooks are the worked examples now, and they are read rather than
/// run: the book links to rendered HTML, and `pages.yml` renders with
/// `nbconvert --to html` and deliberately does not execute. A notebook
/// committed with its outputs cleared would therefore publish as a page of code
/// and no results, and nothing else would notice.
///
/// `nbconvert --execute` stamps every code cell with an integer
/// `execution_count`; a cell that has not been run carries `null`. That is the
/// whole check, and it needs no JSON parser and no interpreter, so it is
/// unconditional.
#[test]
fn every_notebook_was_executed() {
    let dir = repo("notebooks");
    let mut found = 0;
    for entry in fs::read_dir(&dir).expect("notebooks/") {
        let path = entry.expect("a directory entry").path();
        if path.extension().is_none_or(|e| e != "ipynb") {
            continue;
        }
        let name = path
            .file_name()
            .expect("a name")
            .to_string_lossy()
            .into_owned();
        let source = fs::read_to_string(&path).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert!(
            source.contains("\"nbformat\""),
            "notebooks/{name} has no nbformat field, so it is not a notebook"
        );
        let code_cells = source.matches("\"cell_type\": \"code\"").count();
        assert!(
            code_cells > 0,
            "notebooks/{name} has no code cells, which cannot be right"
        );
        assert_eq!(
            source.matches("\"execution_count\": null").count(),
            0,
            "notebooks/{name} has code cells that were never run, so it would \
             publish as code with no results"
        );
        found += 1;
    }
    assert!(found >= 4, "found only {found} notebooks under notebooks/");
}

/// And every notebook still executes against the current build.
///
/// The expensive half, and the one that catches an API change the notebooks
/// have not caught up with. It needs the wheel *and* `jupyter`, and
/// `cargo xtask python` installs only the wheel, so a machine that has not run
/// `pip install jupyter matplotlib` into that virtualenv skips with a message
/// saying so rather than failing.
///
/// Two details are load-bearing, and both were found by this test failing.
///
/// The virtualenv's `bin` goes on `PATH`. `ipykernel` writes a kernelspec whose
/// `argv` starts with the bare word `python`, so the kernel is whatever `PATH`
/// resolves that to, and under `cargo test` that is not the interpreter running
/// nbconvert. The first cell then fails with `No module named 'tepsim'`, which
/// looks exactly like a broken wheel and is not.
///
/// And the working directory is `notebooks/`, not the repository root, because
/// `02-fault-detection-pca.ipynb` does `sys.path.insert(0, str(Path.cwd()))` to
/// import `pcamon` from beside itself.
///
/// `--stdout`, never `--inplace`: re-running the notebooks rewrites every
/// committed output and every plot, which is an author's decision and not a
/// test's.
#[test]
fn every_notebook_still_executes() {
    if !expensive_checks_requested() {
        println!("skipped: set TEP_BOOK=1 to re-execute the notebooks");
        return;
    }
    let Some(python) = venv_python() else {
        println!("skipped: {VENV_PYTHON} does not exist. `cargo xtask python` builds it.");
        return;
    };
    let bin = python.parent().expect("the interpreter has a directory");
    let path = match std::env::var_os("PATH") {
        Some(existing) => {
            let mut dirs = vec![bin.to_path_buf()];
            dirs.extend(std::env::split_paths(&existing));
            std::env::join_paths(dirs).expect("a joinable PATH")
        }
        None => bin.as_os_str().to_os_string(),
    };

    let jupyter = Command::new(&python)
        .args(["-c", "import nbconvert, matplotlib"])
        .current_dir(repo("."))
        .output()
        .expect("the interpreter runs");
    if !jupyter.status.success() {
        println!(
            "skipped: {VENV_PYTHON} has no nbconvert or no matplotlib. \
             `{VENV_PYTHON} -m pip install jupyter matplotlib` adds them; note \
             that `cargo xtask python` rebuilds that virtualenv and drops them."
        );
        return;
    }

    let mut names: Vec<String> = fs::read_dir(repo("notebooks"))
        .expect("notebooks/")
        .filter_map(|e| {
            let path = e.expect("a directory entry").path();
            (path.extension()? == "ipynb").then(|| path.file_name()?.to_str().map(String::from))?
        })
        .collect();
    names.sort();
    let mut ran = 0;
    for name in &names {
        let output = Command::new(&python)
            .args([
                "-m",
                "jupyter",
                "nbconvert",
                "--to",
                "notebook",
                "--execute",
                "--stdout",
                name,
            ])
            .env("PATH", &path)
            .current_dir(repo("notebooks"))
            .output()
            .expect("nbconvert runs");
        assert!(
            output.status.success(),
            "notebooks/{name} no longer executes:\n{}",
            String::from_utf8_lossy(&output.stderr)
        );
        ran += 1;
    }
    assert!(ran >= 4, "executed only {ran} of {names:?}");
}
