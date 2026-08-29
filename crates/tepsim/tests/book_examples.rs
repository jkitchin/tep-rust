//! The book's code listings are the example files, byte for byte.
//!
//! B-0068. `book/src/tutorials/scheduling-a-fault.md` used to carry a
//! `rust,ignore` block and a transcript pasted beside it, and both went stale
//! the moment a default changed: the transcript showed `active [4, 12]` from
//! the driver's forced `IDV(12)`, which stopped happening. Nothing could have
//! caught that, because nothing compiled the listing or ran it.
//!
//! Now the listing is `examples/scheduling_a_fault.rs`, which `cargo test`
//! compiles and this file pins to the page.

use std::fs;
use std::path::{Path, PathBuf};

fn repo(rest: &str) -> PathBuf {
    Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../..")
        .join(rest)
}

/// The body of the first ```rust fence on a page.
fn first_rust_block(markdown: &str) -> String {
    let mut lines = markdown.lines();
    let mut body = Vec::new();
    for line in lines.by_ref() {
        if line.starts_with("```rust") {
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

/// A page, the crate its example lives in, the example's name, and whether the
/// page also quotes that example's output.
///
/// `a-detector` lives in `tepsim-stats` rather than here because it uses that
/// crate, and `xtask ci`'s isolation check forbids `tepsim/Cargo.toml` from so
/// much as naming a development-only crate, dev-dependency included.
///
/// `README.md` shows a snippet and no transcript, so it is pinned for source
/// and skipped for output.
const TUTORIALS: &[(&str, &str, &str, bool)] = &[
    (
        "book/src/tutorials/scheduling-a-fault.md",
        "tepsim",
        "scheduling_a_fault",
        true,
    ),
    (
        "book/src/tutorials/injecting-a-fault.md",
        "tepsim",
        "injecting_a_fault",
        true,
    ),
    (
        "book/src/tutorials/a-detector.md",
        "tepsim-stats",
        "a_detector",
        true,
    ),
    ("README.md", "tepsim", "readme_snippet", false),
];

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

#[test]
fn every_tutorial_shows_the_example_it_runs() {
    for (page, krate, name, _) in TUTORIALS {
        let markdown = fs::read_to_string(repo(page)).unwrap_or_else(|_| panic!("{page}"));
        let example = fs::read_to_string(example_path(krate, name))
            .unwrap_or_else(|_| panic!("{krate}/examples/{name}.rs"));

        // The page omits the file's own `//!` header, which is about the
        // pinning rather than about the tutorial and would be noise on a page.
        let body: String = example
            .lines()
            .skip_while(|l| l.starts_with("//!") || l.is_empty())
            .collect::<Vec<_>>()
            .join("\n");

        let shown = first_rust_block(&markdown);
        // A page that shows no `fn main` is showing its body; see `main_body`.
        let body = if shown.contains("fn main") {
            body
        } else {
            main_body(&body)
        };
        assert!(
            shown.lines().count() >= 4,
            "{page}: the parse found {} lines, so it is reading the wrong fence",
            shown.lines().count()
        );
        assert_eq!(
            shown.trim_end(),
            body.trim_end(),
            "{page} and crates/{krate}/examples/{name}.rs have drifted"
        );
    }
}

/// And the transcript each page quotes is the output that example produces.
///
/// Run rather than trusted. Without this the code could be pinned and the
/// numbers beside it still be from a build two defaults ago, which is exactly
/// what happened: `scheduling-a-fault.md` showed `active [4, 12]` from the
/// driver's forced `IDV(12)` for a day after that stopped being the default.
/// Whether to re-run the examples, which costs about a minute in a debug build.
///
/// Off by default and on in `cargo xtask validate`, in the same spirit as
/// `xtask fidelity` versus the live oracle diff: the cheap check that catches
/// gross breakage runs on every commit, and the expensive one runs
/// periodically. `every_tutorial_shows_the_example_it_runs` is the cheap half
/// and is unconditional, because a listing drifting from its file is the
/// failure that actually happened.
fn transcripts_requested() -> bool {
    std::env::var("TEP_BOOK").is_ok_and(|v| v != "0")
}

#[test]
fn every_tutorial_quotes_the_output_it_gets() {
    if !transcripts_requested() {
        println!(
            "skipped: set TEP_BOOK=1 to re-run the {} tutorial examples \
             (about a minute in a debug build); `cargo xtask validate` does",
            TUTORIALS.len()
        );
        return;
    }
    for (page, krate, name, has_transcript) in TUTORIALS {
        if !has_transcript {
            continue;
        }
        let markdown = fs::read_to_string(repo(page)).unwrap_or_else(|_| panic!("{page}"));
        let quoted = markdown
            .split_once("```text\n")
            .expect("a text fence")
            .1
            .split_once("```")
            .expect("a closing fence")
            .0;

        let output = std::process::Command::new(env!("CARGO"))
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
            .expect("the example runs");
        assert!(
            output.status.success(),
            "{name} failed: {}",
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
            "{page} quotes output {name} no longer produces"
        );
    }
}
