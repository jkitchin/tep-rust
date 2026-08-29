//! The validation figures, drawn as SVG from the runs' own output.
//!
//! [`crate::report`] turns a test binary's stdout into tables. This turns the
//! same stdout into pictures, under the same three rules: only what ran, from
//! the run's own output, and loud rather than clever when the data is not
//! there.
//!
//! # Why the drawing is done here rather than by a library
//!
//! A plotting crate would be the first dependency in this repository that
//! exists for the book rather than for the simulator, and `cargo deny` runs in
//! the gate. The book is built by `mdbook build` alone. A scatter, a band and
//! a log axis are a few hundred lines of `format!`, and this repository
//! already hand-authors SVG in `book/src/images/tep-flowsheet.svg`.
//!
//! # Why the SVG is inlined into the page rather than linked
//!
//! An `<img>` cannot see the host page's colours, so a figure linked that way
//! is legible in exactly one of mdBook's themes and unreadable in the others.
//! Every figure here is included with `{{#include}}`, draws with
//! `currentColor`, and picks its two accents from mdBook's own theme
//! variables: `--links` for the second series and `--warning-border` for a
//! gate. Both are defined in every shipped theme, so a figure that passes in
//! light passes in dark by construction. The flowsheet solved this first.
//!
//! # What a figure is allowed to say
//!
//! Nothing a figure draws is chosen by hand. A dot is orange because its value
//! is on the wrong side of the gate, never because a test's name was
//! recognised, so a genuine regression and a deliberate positive control are
//! drawn the same way and the caption is what tells them apart. That is the
//! property that makes the picture evidence rather than decoration.

// Layout arithmetic, not model arithmetic. `origin + index * pitch` is the
// shape of every row placement here, and rewriting it as `mul_add` would fuse
// a rounding step nobody can see at a tenth of a pixel while making the
// geometry unreadable. The figures are committed to the repository, so the
// expressions also have to keep meaning exactly what they say from one run to
// the next.
#![allow(
    clippy::suboptimal_flops,
    reason = "pixel coordinates: fusing them buys no accuracy and costs clarity"
)]

use std::fmt::Write as _;
use std::path::Path;

use crate::report::{Block, Libm, TargetRun, blocks};

/// Where the figures are written, relative to `book/src`.
pub(crate) const DIR: &str = "book/src/validation/figures";

// ---------------------------------------------------------------------------
// reading numbers out of a transcript
// ---------------------------------------------------------------------------

/// Is this character part of a floating-point literal?
fn numeric(c: char) -> bool {
    c.is_ascii_digit() || matches!(c, '.' | 'e' | 'E' | '+' | '-')
}

/// The first token in `text` that parses as an `f64`.
///
/// Tokens are runs of characters a float can be made of, so `Some(0.075)`
/// yields `0.075` and `at grid#704` yields `704`. Callers slice off everything
/// before the number they want first; see [`float_after`].
pub(crate) fn first_float(text: &str) -> Option<f64> {
    text.split(|c: char| !numeric(c))
        .find_map(|token| token.parse::<f64>().ok())
}

/// The first number after `marker`, or `None` if the marker is absent.
pub(crate) fn float_after(line: &str, marker: &str) -> Option<f64> {
    first_float(line.split_once(marker)?.1)
}

/// `0:9987490 1:12 >=16:7950` as ordered `(bucket, count)` pairs.
///
/// The bucket key is kept as text because the last one is `>=16` rather than a
/// number, and rewriting it as `16` would claim a precision the histogram does
/// not have.
pub(crate) fn ulp_buckets(value: &str) -> Vec<(String, u64)> {
    value
        .split_whitespace()
        .filter_map(|pair| {
            let (bucket, count) = pair.rsplit_once(':')?;
            Some((bucket.to_string(), count.parse::<u64>().ok()?))
        })
        .collect()
}

/// Sort key for a ULP bucket: the number in it, with `>=n` after `n`.
fn bucket_order(bucket: &str) -> (u64, u8) {
    let open = bucket.starts_with(">=");
    let n = first_float(bucket).unwrap_or(0.0).max(0.0) as u64;
    (n, u8::from(open))
}

// ---------------------------------------------------------------------------
// scales
// ---------------------------------------------------------------------------

/// A base-ten logarithmic axis, in pixels.
///
/// Zero has no logarithm, and a validation figure in this project is mostly
/// zeros: every Tier 1 sweep and every `libm-system` comparison is exactly
/// bit-identical. Dropping those points would delete the result, and drawing
/// them at the axis minimum would claim they were merely small. So a zero is
/// drawn in its own lane, off the end of the axis and behind a break, by
/// [`strip`]; the scale itself only handles positive values.
#[derive(Clone, Copy, Debug)]
pub(crate) struct LogAxis {
    /// Exponent at `x0`.
    pub(crate) lo: f64,
    /// Exponent at `x1`.
    pub(crate) hi: f64,
    /// Pixel of `lo`.
    pub(crate) x0: f64,
    /// Pixel of `hi`.
    pub(crate) x1: f64,
}

impl LogAxis {
    /// Where a positive value falls, clamped to the axis.
    pub(crate) fn at(&self, value: f64) -> f64 {
        let exponent = if value > 0.0 {
            value.log10().clamp(self.lo, self.hi)
        } else {
            self.lo
        };
        let fraction = (exponent - self.lo) / (self.hi - self.lo);
        self.x0 + fraction * (self.x1 - self.x0)
    }

    /// Whole-decade tick exponents, thinned so labels never collide.
    pub(crate) fn ticks(&self) -> Vec<f64> {
        let decades = (self.hi - self.lo).max(1.0);
        // About 44 px per label is the narrowest `1e-12` reads at 10 px. The
        // absolute value matters: a vertical axis runs from a large pixel to a
        // small one, and without it every such axis silently got two ticks.
        let wanted = ((self.x1 - self.x0).abs() / 44.0).max(2.0);
        let every = (decades / wanted).ceil().max(1.0) as i64;
        let mut out = Vec::new();
        let mut e = self.lo.ceil() as i64;
        while (e as f64) <= self.hi {
            if e.rem_euclid(every) == 0 {
                out.push(e as f64);
            }
            e += 1;
        }
        out
    }

    /// An axis covering `values` and `must_include`, on whole decades.
    pub(crate) fn covering(values: &[f64], must_include: &[f64], x0: f64, x1: f64) -> Self {
        let mut lo = f64::INFINITY;
        let mut hi = f64::NEG_INFINITY;
        for value in values.iter().chain(must_include).filter(|v| **v > 0.0) {
            let e = value.log10();
            lo = lo.min(e);
            hi = hi.max(e);
        }
        if !lo.is_finite() || !hi.is_finite() {
            // Nothing positive to plot. Six decades ending at one, so the
            // figure still has an axis to hang a zero lane and a gate off.
            lo = -6.0;
            hi = 0.0;
        }
        let mut lo = lo.floor() - 1.0;
        let mut hi = hi.ceil() + 1.0;
        if hi - lo < 3.0 {
            lo -= 1.0;
            hi += 1.0;
        }
        Self { lo, hi, x0, x1 }
    }
}

/// `1e-12`, `1`, `1e3`: an exponent as an axis label.
fn decade_label(exponent: f64) -> String {
    let e = exponent as i64;
    match e {
        0 => "1".to_string(),
        1 => "10".to_string(),
        _ => format!("1e{e}"),
    }
}

/// A measured value, as a figure prints it in a tooltip.
fn value_label(value: f64) -> String {
    if value == 0.0 {
        "exactly 0".to_string()
    } else {
        format!("{value:.3e}")
    }
}

// ---------------------------------------------------------------------------
// the canvas
// ---------------------------------------------------------------------------

/// One figure under construction.
pub(crate) struct Svg {
    id: String,
    width: f64,
    height: f64,
    title: String,
    desc: String,
    body: String,
}

/// `&`, `<` and `>` inside a text node or an attribute.
fn xml(text: &str) -> String {
    text.replace('&', "&amp;")
        .replace('<', "&lt;")
        .replace('>', "&gt;")
        .replace('"', "&quot;")
}

impl Svg {
    fn new(id: &str, width: f64, height: f64, title: &str, desc: &str) -> Self {
        Self {
            id: id.to_string(),
            width,
            height,
            title: title.to_string(),
            desc: desc.to_string(),
            body: String::new(),
        }
    }

    fn line(&mut self, class: &str, x1: f64, y1: f64, x2: f64, y2: f64) {
        let _ = writeln!(
            self.body,
            r#"<line class="{class}" x1="{x1:.1}" y1="{y1:.1}" x2="{x2:.1}" y2="{y2:.1}"/>"#
        );
    }

    fn rect(&mut self, class: &str, x: f64, y: f64, w: f64, h: f64) {
        let _ = writeln!(
            self.body,
            r#"<rect class="{class}" x="{x:.1}" y="{y:.1}" width="{:.1}" height="{h:.1}"/>"#,
            w.max(0.0)
        );
    }

    fn text(&mut self, class: &str, x: f64, y: f64, anchor: &str, body: &str) {
        let anchor = if anchor == "start" {
            String::new()
        } else {
            format!(r#" text-anchor="{anchor}""#)
        };
        let _ = writeln!(
            self.body,
            r#"<text class="{class}" x="{x:.1}" y="{y:.1}"{anchor}>{}</text>"#,
            xml(body)
        );
    }

    /// A data point, with the tooltip that says what it is.
    fn dot(&mut self, class: &str, x: f64, y: f64, r: f64, tip: &str) {
        let _ = writeln!(
            self.body,
            r#"<circle class="{class}" cx="{x:.1}" cy="{y:.1}" r="{r:.1}"><title>{}</title></circle>"#,
            xml(tip)
        );
    }

    /// The whole document, with the provenance the index reads back.
    ///
    /// `provenance` goes in an SVG `<metadata>` element on the second line.
    /// See [`write_figure`] for why it is not the HTML comment the generated
    /// markdown pages use, and why it is not a processing instruction either.
    fn finish(&self, provenance: &str) -> String {
        let id = &self.id;
        format!(
            "<svg id=\"{id}\" viewBox=\"0 0 {:.0} {:.0}\" role=\"img\" \
             aria-labelledby=\"{id}-t {id}-d\" xmlns=\"http://www.w3.org/2000/svg\">\n\
             <metadata>{}</metadata>\n\
             <title id=\"{id}-t\">{}</title><desc id=\"{id}-d\">{}</desc><style>\n{}\
             </style>\n{}</svg>\n",
            self.width,
            self.height,
            xml(provenance),
            xml(&self.title),
            xml(&self.desc),
            css(id),
            self.body
        )
    }
}

/// The stylesheet, scoped to one figure's id.
///
/// Three colours and no more. `currentColor` is the page's own text colour, so
/// it inverts with the theme for free. `--links` and `--warning-border` are
/// mdBook's, defined in light, rust, coal, navy and ayu alike, with fallbacks
/// for anything that renders the file outside the book.
fn css(id: &str) -> String {
    format!(
        "#{id} {{ display: block; width: 100%; height: auto; max-width: 100%; \
         margin: 0 auto; font-family: inherit; }}\n\
         #{id} text {{ fill: currentColor; }}\n\
         #{id} .hd {{ font-size: 13px; font-weight: 600; }}\n\
         #{id} .sub {{ font-size: 11px; fill-opacity: 0.7; }}\n\
         #{id} .lb {{ font-size: 11px; }}\n\
         #{id} .tick {{ font-size: 10px; fill-opacity: 0.72; }}\n\
         #{id} .note {{ font-size: 10.5px; fill-opacity: 0.66; font-style: italic; }}\n\
         #{id} .ax {{ stroke: currentColor; stroke-opacity: 0.5; stroke-width: 1; }}\n\
         #{id} .gr {{ stroke: currentColor; stroke-opacity: 0.16; stroke-width: 1; }}\n\
         #{id} .brk {{ stroke: currentColor; stroke-opacity: 0.38; stroke-width: 1; \
         stroke-dasharray: 2 3; }}\n\
         #{id} .zone {{ fill: currentColor; fill-opacity: 0.055; }}\n\
         #{id} .ok {{ fill: currentColor; fill-opacity: 0.78; }}\n\
         #{id} .okb {{ fill: currentColor; fill-opacity: 0.62; }}\n\
         #{id} .alt {{ fill: var(--links, #20609f); fill-opacity: 0.95; }}\n\
         #{id} .altr {{ fill: none; stroke: var(--links, #20609f); stroke-width: 1.7; }}\n\
         #{id} .bad {{ fill: var(--warning-border, #ff8e00); }}\n\
         #{id} .gate {{ stroke: var(--warning-border, #ff8e00); stroke-width: 1.6; \
         stroke-dasharray: 5 4; }}\n\
         #{id} .band {{ fill: var(--warning-border, #ff8e00); fill-opacity: 0.07; }}\n\
         #{id} .gatetx {{ font-size: 10.5px; fill: var(--warning-border, #ff8e00); \
         font-weight: 600; }}\n"
    )
}

// ---------------------------------------------------------------------------
// figure 1: every comparison against its gate
// ---------------------------------------------------------------------------

/// One comparison block, reduced to what the strip plot draws.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct Point {
    /// The test target the block came from: one lane per target.
    pub(crate) lane: String,
    /// The block's own label, for the tooltip.
    pub(crate) label: String,
    /// The test function, for the tooltip.
    pub(crate) test: String,
    /// Which `libm` the port was built against.
    pub(crate) libm: Libm,
    /// The measured maximum relative error.
    pub(crate) value: f64,
}

/// The comparisons that fell outside the gate, named, without duplicates.
///
/// Naming them in the caption rather than describing them in prose is what
/// makes the figure checkable. There are always a few: a tier reports some
/// comparisons it deliberately does not gate, and Tier 1 injects a
/// mis-transcribed constant on purpose. If a run ever produces a different
/// list, the caption changes with it, which is exactly the signal a sentence
/// written once by hand would have swallowed.
pub(crate) fn outside_the_gate(points: &[Point], gate: f64) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for point in points.iter().filter(|p| p.value > gate) {
        if !out.contains(&point.label) {
            out.push(point.label.clone());
        }
    }
    out
}

/// Every `max rel err` a set of runs printed.
pub(crate) fn error_points(runs: &[TargetRun]) -> Vec<Point> {
    let mut out = Vec::new();
    for run in runs {
        for block in blocks(&run.transcript) {
            let Some(value) = field(&block, "max rel err").and_then(first_float) else {
                continue;
            };
            out.push(Point {
                lane: run.target.clone(),
                label: block.label.clone(),
                test: block.test.clone(),
                libm: run.libm,
                value,
            });
        }
    }
    out
}

/// One field of a block, by key.
fn field<'a>(block: &'a Block, key: &str) -> Option<&'a str> {
    block
        .fields
        .iter()
        .find(|(k, _)| k == key)
        .map(|(_, v)| v.as_str())
}

/// Distinct lanes, in the order the runs produced them.
fn lanes(points: &[Point]) -> Vec<String> {
    let mut out: Vec<String> = Vec::new();
    for point in points {
        if !out.contains(&point.lane) {
            out.push(point.lane.clone());
        }
    }
    out
}

/// The strip plot: every comparison in a tier, against the tier's gate.
///
/// One lane per test target, one dot per comparison, a logarithmic axis, and a
/// lane off the left end for the comparisons that are exactly zero. The gate
/// is a dashed line with everything beyond it shaded, so "did anything cross
/// it" is answered by looking rather than by reading fourteen decimal places.
/// A tier every one of whose comparisons is exactly zero gets no strip plot at
/// all, and [`crate::tier_figures`] says so on the page rather than drawing a
/// logarithmic axis with nothing on it. Tier 3 is that tier: six comparisons,
/// all bit-identical. The ULP figure states the same result in the units that
/// suit it, and an empty error axis with a large shaded failure region reads
/// as a warning about data that does not exist.
pub(crate) fn strip(tier: u8, gate: f64, points: &[Point]) -> Option<Svg> {
    if points.is_empty() || points.iter().all(|p| p.value == 0.0) {
        return None;
    }
    let id = format!("tier{tier}-errors");
    let lanes = lanes(points);
    let values: Vec<f64> = points.iter().map(|p| p.value).collect();

    const WIDTH: f64 = 880.0;
    const GUTTER: f64 = 168.0;
    const ZERO_X: f64 = 196.0;
    const LANE_H: f64 = 32.0;
    let axis = LogAxis::covering(&values, &[gate], 246.0, WIDTH - 24.0);

    // Four text rows above the plot: title, subtitle, the axis and gate
    // labels, then the tick row. They were three, and the gate label ran into
    // the subtitle on every tier.
    let top = 88.0;
    let bottom = top + lanes.len() as f64 * LANE_H;
    let height = bottom + 56.0;

    let zeros = points.iter().filter(|p| p.value == 0.0).count();
    let over = points.iter().filter(|p| p.value > gate).count();
    let platform = points.iter().any(|p| p.libm == Libm::Platform);

    let mut svg = Svg::new(
        &id,
        WIDTH,
        height,
        &format!("Tier {tier}: every comparison against the {gate:e} gate"),
        &format!(
            "A logarithmic strip plot. Each dot is one comparison against the \
             Fortran, placed at its maximum relative error, in a lane named \
             for the test target that produced it. {} of {} comparisons are \
             exactly zero and sit in the separate lane at the left. The gate \
             at {gate:e} is a dashed vertical line with the region beyond it \
             shaded; {over} dots lie beyond it.",
            zeros,
            points.len()
        ),
    );

    svg.text(
        "hd",
        16.0,
        20.0,
        "start",
        &format!("Tier {tier}: every comparison, against the gate"),
    );
    svg.text(
        "sub",
        16.0,
        38.0,
        "start",
        &format!(
            "{} comparisons over {} target(s). {zeros} are bit-identical to \
             the Fortran. {}",
            points.len(),
            lanes.len(),
            match over {
                0 => "None lies beyond the gate.".to_string(),
                1 => "1 lies beyond the gate.".to_string(),
                n => format!("{n} lie beyond the gate."),
            }
        ),
    );

    // The axis, at the top, so a tall figure still shows its scale.
    for exponent in axis.ticks() {
        let x = axis.at(10f64.powf(exponent));
        svg.line("gr", x, top - 8.0, x, bottom);
        svg.text("tick", x, top - 14.0, "middle", &decade_label(exponent));
    }
    svg.line("ax", axis.x0 - 6.0, top - 8.0, axis.x1, top - 8.0);
    svg.text(
        "tick",
        axis.x1,
        top - 30.0,
        "end",
        "maximum relative error against the Fortran",
    );

    // The zero lane, behind a break, because zero is not a small number and
    // drawing it at the axis minimum would say it was merely a small one.
    svg.rect("zone", ZERO_X - 26.0, top - 8.0, 52.0, bottom - top + 8.0);
    svg.text("tick", ZERO_X, top - 14.0, "middle", "= 0");
    svg.line("brk", ZERO_X + 36.0, top - 8.0, ZERO_X + 36.0, bottom);

    // The gate, and everything past it.
    let gate_x = axis.at(gate);
    svg.rect(
        "band",
        gate_x,
        top - 8.0,
        axis.x1 - gate_x,
        bottom - top + 8.0,
    );
    svg.line("gate", gate_x, top - 8.0, gate_x, bottom);
    svg.text(
        "gatetx",
        gate_x - 5.0,
        top - 30.0,
        "end",
        &format!("gate {gate:e}"),
    );

    for (index, lane) in lanes.iter().enumerate() {
        let y = top + index as f64 * LANE_H + LANE_H / 2.0;
        if index % 2 == 1 {
            svg.rect("zone", 8.0, y - LANE_H / 2.0, WIDTH - 16.0, LANE_H);
        }
        svg.text("lb", GUTTER, y + 4.0, "end", lane);
        let mut k = 0usize;
        for point in points.iter().filter(|p| &p.lane == lane) {
            // Deterministic vertical spread, so overlapping dots stay
            // countable. Nothing about a dot's height means anything.
            let offset = f64::from(i32::try_from(k % 3).unwrap_or(0) - 1) * 7.0;
            k += 1;
            let x = if point.value == 0.0 {
                ZERO_X + f64::from(i32::try_from(k % 3).unwrap_or(0) - 1) * 7.0
            } else {
                axis.at(point.value)
            };
            let class = if point.value > gate {
                "bad"
            } else if point.libm == Libm::Platform {
                "alt"
            } else {
                "ok"
            };
            svg.dot(
                class,
                x,
                y + offset,
                3.5,
                &format!(
                    "{}: {} ({} libm, {})",
                    point.label,
                    value_label(point.value),
                    point.libm.label(),
                    point.test
                ),
            );
        }
    }

    // Legend.
    let y = bottom + 26.0;
    let mut x = 16.0;
    svg.dot("ok", x, y - 4.0, 3.5, "vendored libm");
    svg.text("note", x + 9.0, y, "start", "vendored libm");
    x += 130.0;
    if platform {
        svg.dot("alt", x, y - 4.0, 3.5, "platform libm");
        svg.text(
            "note",
            x + 9.0,
            y,
            "start",
            "platform libm, bit-exact claim",
        );
        x += 210.0;
    }
    svg.dot("bad", x, y - 4.0, 3.5, "beyond the gate");
    svg.text("note", x + 9.0, y, "start", "beyond the gate");
    svg.text(
        "note",
        WIDTH - 16.0,
        y,
        "end",
        "hover a dot for the comparison it came from",
    );
    Some(svg)
}

// ---------------------------------------------------------------------------
// figure 2: how many bits actually differ
// ---------------------------------------------------------------------------

/// The ULP histogram of a tier, aggregated per `libm`.
///
/// Blocks beyond the gate are left out and the figure says so: they are the
/// positive controls, whose whole purpose is to be enormous, and one of them
/// contributes ten million comparisons at `>=16` ULP that would otherwise be
/// the only thing this figure showed.
pub(crate) fn ulp(tier: u8, gate: f64, runs: &[TargetRun]) -> Option<Svg> {
    let mut per_libm: Vec<(Libm, Vec<(String, u64)>)> = Vec::new();
    let mut excluded = 0usize;
    for run in runs {
        for block in blocks(&run.transcript) {
            let Some(histogram) = field(&block, "ulp histogram") else {
                continue;
            };
            let error = field(&block, "max rel err").and_then(first_float);
            if error.is_some_and(|e| e > gate) {
                excluded += 1;
                continue;
            }
            let slot = match per_libm.iter_mut().find(|(l, _)| *l == run.libm) {
                Some((_, counts)) => counts,
                None => {
                    per_libm.push((run.libm, Vec::new()));
                    &mut per_libm.last_mut()?.1
                }
            };
            for (bucket, count) in ulp_buckets(histogram) {
                match slot.iter_mut().find(|(b, _)| *b == bucket) {
                    Some((_, total)) => *total += count,
                    None => slot.push((bucket, count)),
                }
            }
        }
    }
    per_libm.retain(|(_, counts)| !counts.is_empty());
    if per_libm.is_empty() {
        return None;
    }
    for (_, counts) in &mut per_libm {
        counts.sort_by_key(|(bucket, _)| bucket_order(bucket));
    }

    let id = format!("tier{tier}-ulp");
    const WIDTH: f64 = 880.0;
    const ROW_H: f64 = 20.0;
    let rows: usize = per_libm.iter().map(|(_, c)| c.len()).sum();
    let top = 72.0;
    // Each group is a header, its bars, and then its own tick row: 26 above
    // and 26 below. Counting only the header put the last group's tick labels
    // on top of the closing note.
    let height = top + rows as f64 * ROW_H + per_libm.len() as f64 * 52.0 + 34.0;

    let biggest = per_libm
        .iter()
        .flat_map(|(_, c)| c.iter().map(|(_, n)| *n as f64))
        .fold(1.0_f64, f64::max);
    // Built rather than fitted: a count axis starts at one, because a bar
    // running from a tenth of a comparison would be a made-up quantity, and
    // the bars are drawn from `x0` on the assumption that `x0` is one.
    let axis = LogAxis {
        lo: 0.0,
        hi: biggest.log10().ceil().max(1.0) + 0.25,
        x0: 176.0,
        x1: WIDTH - 96.0,
    };
    let total: u64 = per_libm
        .iter()
        .flat_map(|(_, c)| c.iter().map(|(_, n)| *n))
        .sum();

    let mut svg = Svg::new(
        &id,
        WIDTH,
        height,
        &format!("Tier {tier}: how many bits differ"),
        &format!(
            "A bar chart on a logarithmic count axis. Each bar is the number \
             of comparisons whose result differed from the Fortran's by that \
             many units in the last place, grouped by which libm the port was \
             built against. {} comparisons in total.",
            group(total)
        ),
    );
    svg.text(
        "hd",
        16.0,
        20.0,
        "start",
        &format!("Tier {tier}: how many bits actually differ"),
    );
    svg.text(
        "sub",
        16.0,
        36.0,
        "start",
        &format!(
            "{} comparisons, by units in the last place. The count axis is \
             logarithmic.{}",
            group(total),
            match excluded {
                0 => String::new(),
                1 => " One block beyond the gate is left out.".to_string(),
                n => format!(" {n} blocks beyond the gate are left out."),
            }
        ),
    );

    let mut y = top;
    for (libm, counts) in &per_libm {
        let subtotal: u64 = counts.iter().map(|(_, n)| *n).sum();
        let zero = counts.iter().find(|(b, _)| b == "0").map_or(0, |(_, n)| *n);
        svg.text(
            "lb",
            16.0,
            y + 12.0,
            "start",
            &format!(
                "{} libm: {} of {} identical to the last bit",
                libm.label(),
                group(zero),
                group(subtotal)
            ),
        );
        y += 26.0;
        for exponent in axis.ticks() {
            let x = axis.at(10f64.powf(exponent));
            svg.line("gr", x, y - 4.0, x, y + counts.len() as f64 * ROW_H - 2.0);
        }
        for (bucket, count) in counts {
            let label = if bucket == "1" {
                "1 ULP".to_string()
            } else {
                format!("{bucket} ULP")
            };
            svg.text("tick", 168.0, y + 13.0, "end", &label);
            let x = axis.at(*count as f64);
            let class = if bucket == "0" { "ok" } else { "okb" };
            svg.rect(class, axis.x0, y + 3.0, x - axis.x0, ROW_H - 8.0);
            svg.text("tick", x + 6.0, y + 13.0, "start", &group(*count));
            y += ROW_H;
        }
        for exponent in axis.ticks() {
            let x = axis.at(10f64.powf(exponent));
            svg.text("tick", x, y + 13.0, "middle", &decade_label(exponent));
        }
        svg.line("ax", axis.x0, y + 1.0, axis.x1, y + 1.0);
        y += 26.0;
    }
    svg.text(
        "note",
        16.0,
        height - 12.0,
        "start",
        "A bar exists only where the count is not zero, so a tier with one bar \
         differed nowhere.",
    );
    Some(svg)
}

/// `9987490` as `9,987,490`.
fn group(n: u64) -> String {
    let digits = n.to_string();
    let mut out = String::new();
    for (index, digit) in digits.chars().enumerate() {
        if index > 0 && (digits.len() - index) % 3 == 0 {
            out.push(',');
        }
        out.push(digit);
    }
    out
}

// ---------------------------------------------------------------------------
// figure 3: Tier 4 against the instrument noise
// ---------------------------------------------------------------------------

/// One scenario's trajectory error, as a fraction of the instrument noise.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct NoisePoint {
    /// `nominal`, `IDV(7)`.
    pub(crate) scenario: String,
    /// Which `libm` produced it.
    pub(crate) libm: Libm,
    /// Worst error over the run, divided by `XNS(i)` for that measurement.
    pub(crate) worst: f64,
    /// How long the two stayed inside the noise band, in hours.
    pub(crate) hours: f64,
}

/// Read the Tier 4 scenario sweep out of a transcript.
///
/// The line shape is the one `tier4_fault_scenarios` prints:
/// `IDV(10)    within XNS 4.000 h, first out none, worst 3.66e-8`.
///
/// All three markers are required, and `first out` is the one that earns its
/// keep. Another test in the same binary prints `vendored libm: within XNS for
/// 4.000 h, worst 2.446e-5 of XNS at 4 h`, which matched on two markers and put
/// a row labelled `test` at the top of the figure, from libtest's own prefix.
pub(crate) fn noise_points(runs: &[TargetRun]) -> Vec<NoisePoint> {
    let mut out = Vec::new();
    for run in runs {
        for line in &run.transcript {
            if !line.contains("within XNS ")
                || !line.contains(", first out ")
                || !line.contains(", worst ")
            {
                continue;
            }
            let trimmed = line.trim();
            let Some(scenario) = trimmed.split_whitespace().next() else {
                continue;
            };
            // The nominal 8 h block prints `within XNS for : 8.000 h` on its
            // own line and its worst on another, so it never matches here.
            let (Some(hours), Some(worst)) = (
                float_after(trimmed, "within XNS "),
                float_after(trimmed, ", worst "),
            ) else {
                continue;
            };
            out.push(NoisePoint {
                scenario: scenario.to_string(),
                libm: run.libm,
                worst,
                hours,
            });
        }
    }
    out
}

/// The sweep's horizon, from the header the sweep prints above its rows.
///
/// `Tier 4, fault scenarios, 4 h each:`. Read rather than assumed, so that
/// changing the test's horizon changes the figure's caption with it.
pub(crate) fn sweep_hours(runs: &[TargetRun]) -> Option<f64> {
    runs.iter()
        .flat_map(|run| run.transcript.iter())
        .find_map(|line| float_after(line, "fault scenarios, "))
}

/// Trajectory error against the band the instruments could resolve.
///
/// The claim Tier 4 makes is not that the two trajectories agree, which after
/// eight simulated hours of different `exp` implementations they cannot. It is
/// that the disagreement is far below what the plant's own instruments could
/// see. So the figure's reference line is `XNS(i)`, the measurement noise
/// standard deviation, and every point is an error divided by it.
pub(crate) fn noise_band(points: &[NoisePoint], hours: f64) -> Option<Svg> {
    if points.is_empty() {
        return None;
    }
    let id = "tier4-noise-band";
    let mut scenarios: Vec<String> = Vec::new();
    for point in points {
        if !scenarios.contains(&point.scenario) {
            scenarios.push(point.scenario.clone());
        }
    }

    const WIDTH: f64 = 880.0;
    const ROW_H: f64 = 19.0;
    const ZERO_X: f64 = 178.0;
    let values: Vec<f64> = points.iter().map(|p| p.worst).collect();
    let axis = LogAxis::covering(&values, &[1.0], 228.0, WIDTH - 130.0);

    // Four text rows above the plot, as in `strip`: title, subtitle, the axis
    // and band labels, then the ticks.
    let top = 88.0;
    let bottom = top + scenarios.len() as f64 * ROW_H;
    let height = bottom + 74.0;

    let zeros = points.iter().filter(|p| p.worst == 0.0).count();
    let outside = points.iter().filter(|p| p.worst >= 1.0).count();
    let worst = values.iter().copied().fold(0.0_f64, f64::max);

    let mut svg = Svg::new(
        id,
        WIDTH,
        height,
        "Tier 4: trajectory error against the instrument noise",
        &format!(
            "A logarithmic strip plot with one row per scenario. Each marker \
             is the worst disagreement between the port and the Fortran over a \
             {hours} hour run, divided by the measurement noise standard \
             deviation of the channel it occurred on. The band at one, where \
             the error would equal the instrument noise, is shaded; {outside} \
             of {} markers reach it. {zeros} are exactly zero.",
            points.len()
        ),
    );
    svg.text(
        "hd",
        16.0,
        20.0,
        "start",
        "Tier 4: how far apart the trajectories get, in units of instrument noise",
    );
    svg.text(
        "sub",
        16.0,
        36.0,
        "start",
        &format!(
            "{} scenarios, {hours} h each. Worst error over the run divided by \
             XNS(i). {zeros} of {} runs are bit-identical.",
            scenarios.len(),
            points.len()
        ),
    );

    for exponent in axis.ticks() {
        let x = axis.at(10f64.powf(exponent));
        svg.line("gr", x, top - 8.0, x, bottom);
        svg.text("tick", x, top - 14.0, "middle", &decade_label(exponent));
    }
    svg.line("ax", axis.x0 - 6.0, top - 8.0, axis.x1, top - 8.0);
    svg.text(
        "tick",
        axis.x0,
        top - 32.0,
        "start",
        "worst error over the run, as a fraction of that channel's noise",
    );

    svg.rect("zone", ZERO_X - 24.0, top - 8.0, 48.0, bottom - top + 8.0);
    svg.text("tick", ZERO_X, top - 14.0, "middle", "= 0");
    svg.line("brk", ZERO_X + 34.0, top - 8.0, ZERO_X + 34.0, bottom);

    // The noise band: at one the error equals the instrument noise. The label
    // hangs to the left of its line, because the line sits near the right edge
    // and a label started there runs off the canvas.
    let band_x = axis.at(1.0);
    svg.rect(
        "band",
        band_x,
        top - 8.0,
        axis.x1 - band_x,
        bottom - top + 8.0,
    );
    svg.line("gate", band_x, top - 8.0, band_x, bottom);
    svg.text("gatetx", band_x - 6.0, top - 32.0, "end", "1 x XNS(i)");

    for (index, scenario) in scenarios.iter().enumerate() {
        let y = top + index as f64 * ROW_H + ROW_H / 2.0;
        if index % 2 == 1 {
            svg.rect("zone", 8.0, y - ROW_H / 2.0, WIDTH - 16.0, ROW_H);
        }
        svg.text("lb", 150.0, y + 4.0, "end", scenario);
        for point in points.iter().filter(|p| &p.scenario == scenario) {
            let x = if point.worst == 0.0 {
                ZERO_X
            } else {
                axis.at(point.worst)
            };
            let tip = format!(
                "{scenario}, {} libm: worst {}, inside the noise for {:.3} h",
                point.libm.label(),
                value_label(point.worst),
                point.hours
            );
            match point.libm {
                Libm::Platform => {
                    let _ = writeln!(
                        svg.body,
                        r#"<circle class="altr" cx="{x:.1}" cy="{y:.1}" r="4.2"><title>{}</title></circle>"#,
                        xml(&tip)
                    );
                }
                Libm::Vendored => svg.dot("ok", x, y, 3.6, &tip),
            }
        }
    }

    let y = bottom + 26.0;
    svg.dot("ok", 16.0, y - 4.0, 3.6, "vendored libm");
    svg.text("note", 25.0, y, "start", "vendored libm");
    let _ = writeln!(
        svg.body,
        r#"<circle class="altr" cx="146.0" cy="{:.1}" r="4.2"/>"#,
        y - 4.0
    );
    svg.text(
        "note",
        156.0,
        y,
        "start",
        "platform libm: the same exp gfortran calls",
    );
    svg.text(
        "note",
        WIDTH - 16.0,
        y,
        "end",
        &format!("worst anywhere: {}", value_label(worst)),
    );
    svg.text(
        "note",
        16.0,
        y + 18.0,
        "start",
        "The shaded band is where the disagreement would be as large as the \
         channel's own measurement noise.",
    );
    Some(svg)
}

// ---------------------------------------------------------------------------
// figure 4: Tier 5, cross-source against the within-source null
// ---------------------------------------------------------------------------

/// One statistic, measured across the two sources and inside the reference.
#[derive(Clone, Debug, PartialEq)]
pub(crate) struct CalibrationPoint {
    /// What was measured, and where.
    pub(crate) what: String,
    /// The statistic between the Fortran and the port.
    pub(crate) cross: f64,
    /// The same statistic between two halves of the Fortran's own runs.
    pub(crate) within: f64,
}

/// Read the Tier 5 calibration out of the battery's transcript.
///
/// Two line shapes, both printed by `tier5_battery`:
///
/// ```text
///   frobenius: cross 6.973021e-12, within-source max 1.090611e1, p 1.0000
///   KS cross 1.0000 vs within Some(0.07500000000000001)
/// ```
///
/// The second comes from the test that shifts one variable by ten standard
/// deviations, and it is kept precisely because it is a real difference
/// measured on the same axes.
pub(crate) fn calibration_points(runs: &[TargetRun]) -> Vec<CalibrationPoint> {
    let mut out = Vec::new();
    let mut scenario = String::new();
    for run in runs {
        for line in &run.transcript {
            let trimmed = line.trim();
            if let Some(rest) = trimmed.strip_prefix("=== ") {
                scenario = rest
                    .split_once(" (")
                    .map_or(rest, |(name, _)| name)
                    .to_string();
            }
            if trimmed.contains("frobenius:") {
                if let (Some(cross), Some(within)) = (
                    float_after(trimmed, "cross "),
                    float_after(trimmed, "within-source max "),
                ) {
                    out.push(CalibrationPoint {
                        what: format!("correlation matrix, {scenario}"),
                        cross,
                        within,
                    });
                }
            } else if trimmed.contains(" cross ") && trimmed.contains(" vs within ") {
                let what = trimmed.split_whitespace().next().unwrap_or("statistic");
                if let (Some(cross), Some(within)) = (
                    float_after(trimmed, " cross "),
                    float_after(trimmed, " vs within "),
                ) {
                    out.push(CalibrationPoint {
                        what: format!("{what}, one variable shifted 10 sd"),
                        cross,
                        within,
                    });
                }
            }
        }
    }
    out
}

/// How big a battery produced these numbers, in the battery's own words.
///
/// `Tier 5 battery: 3 scenarios x 4 seeds x 2 h, 40 samples per run, vendored
/// libm`. The smoke battery and the full one differ by two orders of magnitude
/// in every direction, so a figure that did not say which it came from would
/// be unreadable as evidence.
pub(crate) fn battery_size(runs: &[TargetRun]) -> Option<String> {
    runs.iter()
        .flat_map(|run| run.transcript.iter())
        .find_map(|line| Some(line.split_once("Tier 5 battery: ")?.1.trim().to_string()))
}

/// The cross-source value against the null it is judged by.
///
/// This is the claim Tier 5 actually makes, and it is not "the two sources are
/// close". No absolute scale exists for a Kolmogorov-Smirnov statistic or a
/// Frobenius distance, so the battery builds one: split the Fortran's own runs
/// in half, compute the statistic Fortran against Fortran, and use that spread
/// as the yardstick. A point above the diagonal means the two implementations
/// differ from each other by *less* than the reference differs from itself.
pub(crate) fn calibration(points: &[CalibrationPoint], size: &str) -> Option<Svg> {
    if points.is_empty() {
        return None;
    }
    let id = "tier5-calibration";
    const WIDTH: f64 = 880.0;
    const HEIGHT: f64 = 590.0;
    let plot_top = 76.0;
    let plot_bottom = 476.0;

    let mut values: Vec<f64> = Vec::new();
    for point in points {
        values.push(point.cross);
        values.push(point.within);
    }
    // One scale for both axes, so the diagonal really is the locus of
    // equality. A figure whose 45-degree line was not `y = x` would be a lie
    // told with geometry, so the plot is squared off to the smaller of the two
    // available extents rather than stretched to fill the canvas.
    let wide = LogAxis::covering(&values, &[], 300.0, WIDTH - 40.0);
    let decades = wide.hi - wide.lo;
    let side = (wide.x1 - wide.x0).min(plot_bottom - plot_top);
    let axis = LogAxis {
        x1: wide.x0 + side,
        ..wide
    };
    let vertical = LogAxis {
        lo: axis.lo,
        hi: axis.hi,
        x0: plot_bottom,
        x1: plot_bottom - side,
    };
    let span = side;
    debug_assert!(decades > 0.0, "an axis with no decades cannot be squared");

    let inside = points.iter().filter(|p| p.cross <= p.within).count();
    let mut svg = Svg::new(
        id,
        WIDTH,
        HEIGHT,
        "Tier 5: the cross-source difference against the reference's own spread",
        &format!(
            "A logarithmic scatter plot. The horizontal axis is a statistic \
             measured between the Fortran and the Rust port; the vertical axis \
             is the same statistic measured between two halves of the \
             Fortran's own runs. Both axes share one scale, so the diagonal is \
             the line of equality. {inside} of {} points lie above it, meaning \
             the two implementations differ by less than the reference differs \
             from itself.",
            points.len()
        ),
    );
    svg.text(
        "hd",
        16.0,
        20.0,
        "start",
        "Tier 5: the two sources against the reference's own run-to-run spread",
    );
    svg.text(
        "sub",
        16.0,
        36.0,
        "start",
        &format!(
            "{size}. A point above the diagonal differs across sources by less \
             than the reference differs from itself.",
        ),
    );

    // The passing half-plane, drawn as a polygon above the diagonal.
    let _ = writeln!(
        svg.body,
        r#"<path class="zone" d="M{:.1} {:.1} L{:.1} {:.1} L{:.1} {:.1} Z"/>"#,
        axis.x0,
        vertical.x0,
        axis.x0 + side,
        vertical.x1,
        axis.x0,
        vertical.x1
    );

    for exponent in axis.ticks() {
        let x = axis.at(10f64.powf(exponent));
        svg.line("gr", x, plot_top, x, vertical.x0);
        svg.text(
            "tick",
            x,
            vertical.x0 + 16.0,
            "middle",
            &decade_label(exponent),
        );
    }
    for exponent in vertical.ticks() {
        let y = vertical.at(10f64.powf(exponent));
        svg.line("gr", axis.x0, y, axis.x1, y);
        svg.text(
            "tick",
            axis.x0 - 8.0,
            y + 4.0,
            "end",
            &decade_label(exponent),
        );
    }
    svg.line("ax", axis.x0, vertical.x0, axis.x1, vertical.x0);
    svg.line("ax", axis.x0, plot_top, axis.x0, vertical.x0);

    // The diagonal: equality between the two.
    let reach = side;
    svg.line(
        "gate",
        axis.x0,
        vertical.x0,
        axis.x0 + reach,
        vertical.x0 - reach,
    );
    svg.text(
        "gatetx",
        axis.x0 + reach - 6.0,
        vertical.x0 - reach + 16.0,
        "end",
        "equal: the port is as different as a rerun",
    );

    svg.text(
        "tick",
        axis.x0 + span / 2.0,
        vertical.x0 + 34.0,
        "middle",
        "measured between the Fortran and the port",
    );
    svg.text(
        "note",
        axis.x0 - 18.0,
        plot_top + 14.0,
        "end",
        "measured inside",
    );
    svg.text(
        "note",
        axis.x0 - 18.0,
        plot_top + 27.0,
        "end",
        "the Fortran alone",
    );

    // Labels are stacked rather than placed at each point, because the points
    // that matter cluster: the three correlation-matrix results sit within a
    // pixel of one another and their labels drew on top of each other. The
    // leader line is what keeps a stacked label attached to its point.
    let mut used: Vec<(f64, f64)> = Vec::new();
    let mut unlabelled = 0usize;
    for point in points {
        let x = axis.at(point.cross);
        let y = vertical.at(point.within);
        let pass = point.cross <= point.within;
        let tip = format!(
            "{}: across sources {}, inside the Fortran {}",
            point.what,
            value_label(point.cross),
            value_label(point.within)
        );
        svg.dot(if pass { "ok" } else { "bad" }, x, y, 5.5, &tip);

        // Stacked away from the point rather than through it: the passing
        // cluster sits within a pixel of itself, and stacking downwards laid
        // the second label straight across the dots it was naming.
        let step = if pass { -13.0 } else { 13.0 };
        let mut label_y = y + if pass { -11.0 } else { 19.0 };
        while used
            .iter()
            .any(|(ux, uy)| (uy - label_y).abs() < 12.0 && (ux - x).abs() < 260.0)
        {
            label_y += step;
        }
        // A label that has been pushed off the plot is dropped rather than
        // drawn over the axis. The smoke battery has five points and never
        // reaches this; the full one has twenty-two, and a stack of
        // twenty-two labels would bury the figure it is annotating. Every
        // dot keeps its tooltip either way, and the count below says how many
        // went unlabelled.
        used.push((x, label_y));
        if label_y < plot_top + 10.0 || label_y > vertical.x0 - 6.0 {
            unlabelled += 1;
            continue;
        }
        if (label_y - y).abs() > 14.0 {
            let anchor = if label_y < y { label_y } else { label_y - 4.0 };
            svg.line("brk", x + 4.0, y, x + 9.0, anchor);
        }
        svg.text("note", x + 11.0, label_y, "start", &point.what);
    }

    let y = HEIGHT - 46.0;
    svg.dot("ok", 16.0, y - 4.0, 4.5, "inside the null");
    svg.text(
        "note",
        26.0,
        y,
        "start",
        "above the diagonal: inside the reference's own spread",
    );
    svg.dot("bad", 430.0, y - 4.0, 4.5, "outside the null");
    svg.text(
        "note",
        440.0,
        y,
        "start",
        "below it: a difference the calibration can see",
    );
    svg.text(
        "note",
        16.0,
        HEIGHT - 22.0,
        "start",
        &match unlabelled {
            0 => "The points below the diagonal are the battery's own positive \
                  control: one variable of the reference shifted by ten \
                  standard deviations and compared with the unshifted \
                  reference."
                .to_string(),
            n => format!(
                "The points below the diagonal are the battery's own positive \
                 control: one variable of the reference shifted by ten \
                 standard deviations. {n} point(s) are unlabelled for want of \
                 room; hover for the name."
            ),
        },
    );
    Some(svg)
}

// ---------------------------------------------------------------------------
// writing
// ---------------------------------------------------------------------------

/// Write one figure. Its file name is its own id, so a page can find it.
///
/// The provenance is the first line of the file, in the shape
/// [`crate::report::read_provenance`] parses, so a page that inlines a figure
/// some *other* run drew can still say which run that was. Without it the
/// index would have to claim every figure came from the latest run, which is
/// false the moment someone regenerates one tier and not another.
///
/// # Why it is `<metadata>` and not a comment
///
/// The generated markdown pages carry the same line as an HTML comment, and
/// the first version of this did too. Two things ruled that out. An XML
/// comment may not contain `--`, and every command that writes one of these
/// has `--tiers` in it: mdBook inlines the figure into a page, where an HTML
/// parser tolerates that, but it *also* copies the `.svg` to the site, and
/// opening that file directly serves it as `image/svg+xml`, where a strict
/// parser rejects the whole document. A processing instruction fixes that and
/// breaks something else: mdBook parses the page's HTML to place its heading
/// anchors, and `<?` makes it log `html parse error in validation.md`.
/// `<metadata>` is a real SVG element, valid in both parsers, and it renders
/// nothing.
pub(crate) fn write_figure(root: &Path, svg: &Svg, command: &str) -> Result<(), String> {
    let commit = crate::report::describe_commit(root);
    let provenance = format!(
        "GENERATED by `{command}` from commit `{commit}`. Do not edit by \
         hand: the next run overwrites it."
    );
    crate::report::write_generated(
        root,
        &format!("{DIR}/{}.svg", svg.id),
        &svg.finish(&provenance),
    )
}

/// A figure that a page inlines, with everything the caption has to say.
///
/// `caption` states what would make the picture false. A figure without a
/// stated failure condition is decoration, so this is a required field rather
/// than an option, and `from` names the `LOG.org` iteration whose measurement
/// the picture repeats.
#[derive(Clone, Copy, Debug)]
pub(crate) struct Caption<'a> {
    /// The figure's id, which is also its file name.
    pub(crate) id: &'a str,
    /// The bold lead-in, matching the figure's own `<title>`.
    pub(crate) title: &'a str,
    /// What it shows, and what would falsify it.
    pub(crate) caption: &'a str,
    /// The `LOG.org` iteration that first measured this.
    pub(crate) from: &'a str,
}

/// The markdown that inlines a figure, or `None` when no run has drawn it.
///
/// `prefix` is the path from the including page to [`DIR`]: `figures/` for the
/// generated chapters, `validation/figures/` for the narrative page one level
/// up.
///
/// The `<figure>` wrapper is the convention `book/src/introduction.md` already
/// uses for the flowsheet, and it is load-bearing twice over. Without it the
/// SVG lands inside a paragraph, because its opening tag shares a line with
/// the `<title>` and CommonMark only starts an HTML block for a tag alone on
/// its line. And it puts the caption in a `<figcaption>`, where a screen
/// reader will read it as the figure's caption rather than as the next
/// sentence of the page.
pub(crate) fn include(root: &Path, prefix: &str, figure: &Caption<'_>) -> Option<String> {
    let relative = format!("{DIR}/{}.svg", figure.id);
    let (by, commit) = crate::report::read_provenance(root, &relative)?;
    Some(format!(
        "<figure style=\"margin: 1.5rem 0;\">\n\
         {{{{#include {prefix}{}.svg}}}}\n\
         <figcaption style=\"font-size: 0.9em; opacity: 0.8; margin-top: 0.6em;\">\
         <strong>{}.</strong> {} Drawn by <code>{}</code> at commit \
         <code>{}</code>; the measurement it repeats was first recorded in {} \
         (<code>LOG.org</code>).</figcaption>\n\
         </figure>\n",
        figure.id,
        inline_html(figure.title),
        inline_html(figure.caption),
        xml(&by),
        xml(&commit),
        inline_html(figure.from)
    ))
}

/// A caption's markdown as the inline HTML a `<figcaption>` needs.
///
/// mdBook passes a raw HTML block through without parsing markdown inside it,
/// so the two bits of markup the captions use have to be converted here.
/// Writing the captions as HTML in the first place would make them unreadable
/// where they are defined, which is the only place anyone edits them.
fn inline_html(text: &str) -> String {
    let mut out = String::new();
    for (index, chunk) in text.split('`').enumerate() {
        if index % 2 == 1 {
            let _ = write!(out, "<code>{}</code>", xml(chunk));
        } else {
            // Bold, the other thing a caption uses. Only when the markers
            // pair: a stray `**` stays as text rather than bolding the rest of
            // the caption or, worse, silently vanishing from it.
            let runs: Vec<&str> = chunk.split("**").collect();
            if runs.len() % 2 == 0 {
                out.push_str(&xml(chunk));
                continue;
            }
            for (bold, run) in runs.iter().enumerate() {
                if bold % 2 == 1 {
                    let _ = write!(out, "<strong>{}</strong>", xml(run));
                } else {
                    out.push_str(&xml(run));
                }
            }
        }
    }
    out
}

#[cfg(test)]
mod tests {
    // A figure's whole job here is to distinguish "exactly zero" from "small",
    // so the tests compare against zero exactly. A near-equality assertion
    // would pass on a value the figure is required to draw somewhere else.
    #![allow(
        clippy::float_cmp,
        reason = "exact zero is the property under test, not a nearby value"
    )]

    use super::*;

    fn run(target: &str, libm: Libm, transcript: &str) -> TargetRun {
        TargetRun {
            target: target.to_string(),
            command: "cargo test".to_string(),
            libm,
            transcript: transcript.lines().map(str::to_string).collect(),
            tally: None,
        }
    }

    #[test]
    fn a_float_is_found_after_its_marker() {
        let line = "  frobenius: cross 6.973021e-12, within-source max 1.090611e1, p 1.0000";
        assert_eq!(float_after(line, "cross "), Some(6.973_021e-12));
        assert_eq!(float_after(line, "within-source max "), Some(1.090_611e1));
        assert_eq!(float_after(line, "nothing like this"), None);
    }

    /// The battery prints a `Debug`-formatted `Option`, so the number the
    /// figure needs is inside `Some(...)`.
    #[test]
    fn a_number_inside_a_debug_option_is_still_a_number() {
        let line = "  KS cross 1.0000 vs within Some(0.07500000000000001)";
        assert_eq!(float_after(line, " cross "), Some(1.0));
        assert_eq!(
            float_after(line, " vs within "),
            Some(0.075_000_000_000_000_01)
        );
    }

    #[test]
    fn a_value_with_a_location_yields_the_value() {
        assert_eq!(first_float("1.597e-5 at grid#704 T=21.875"), Some(1.597e-5));
        assert_eq!(first_float("0.000e0 at grid#0 T=0"), Some(0.0));
        assert_eq!(first_float("no numbers here"), None);
    }

    #[test]
    fn a_ulp_histogram_parses_into_ordered_buckets() {
        let mut buckets = ulp_buckets("2:3 0:9700 1:12 >=16:5");
        buckets.sort_by_key(|(b, _)| bucket_order(b));
        assert_eq!(
            buckets,
            vec![
                ("0".to_string(), 9700),
                ("1".to_string(), 12),
                ("2".to_string(), 3),
                (">=16".to_string(), 5),
            ]
        );
    }

    #[test]
    fn a_log_axis_places_its_ends_and_its_middle() {
        let axis = LogAxis {
            lo: -4.0,
            hi: 0.0,
            x0: 0.0,
            x1: 400.0,
        };
        assert!((axis.at(1e-4) - 0.0).abs() < 1e-9);
        assert!((axis.at(1.0) - 400.0).abs() < 1e-9);
        assert!((axis.at(1e-2) - 200.0).abs() < 1e-9);
        // Off the end, and off the bottom: clamped rather than drawn outside.
        assert!((axis.at(1e-9) - 0.0).abs() < 1e-9);
        assert!((axis.at(0.0) - 0.0).abs() < 1e-9);
        assert!((axis.at(1e9) - 400.0).abs() < 1e-9);
    }

    /// The gate has to be on the axis even when nothing measured is near it,
    /// which is the usual case here: every Tier 1 sweep is exactly zero.
    #[test]
    fn the_axis_always_covers_the_gate() {
        let axis = LogAxis::covering(&[0.0, 0.0], &[1e-12], 0.0, 100.0);
        assert!(axis.lo <= -12.0 && axis.hi >= -12.0, "{axis:?}");
    }

    #[test]
    fn ticks_thin_out_rather_than_overlapping() {
        let axis = LogAxis {
            lo: -40.0,
            hi: 0.0,
            x0: 0.0,
            x1: 200.0,
        };
        let ticks = axis.ticks();
        assert!(ticks.len() <= 6, "{ticks:?}");
        assert!(ticks.windows(2).all(|w| w[1] > w[0]));
    }

    #[test]
    fn error_points_come_from_the_blocks_a_run_printed() {
        let runs = [run(
            "tier2_kinetics",
            Libm::Vendored,
            "running 1 test\n\
             test the_kinetics_match ... tier1 kinetics RR\n\
             \x20 cases          : 9700\n\
             \x20 max rel err    : 8.492e-16 at perturbed#618[2]\n\
             ok\n",
        )];
        let points = error_points(&runs);
        assert_eq!(points.len(), 1);
        assert_eq!(points[0].lane, "tier2_kinetics");
        assert_eq!(points[0].label, "tier1 kinetics RR");
        assert_eq!(points[0].test, "the_kinetics_match");
        assert!((points[0].value - 8.492e-16).abs() < 1e-30);
    }

    /// The rule that colours a dot is the value against the gate, never a
    /// test's name. Both figures below are generated from the same input, and
    /// the only orange dot is the one that is actually outside.
    #[test]
    fn a_dot_is_orange_because_of_its_value_and_nothing_else() {
        let runs = [run(
            "tier1_enthalpy",
            Libm::Vendored,
            "test whatever ... tier1 TESUB1 ity=2, offset as f64\n\
             \x20 max rel err    : 1.597e-5 at grid#704\n\
             \x20 ulp histogram  : >=16:7950\n\
             tier1 TESUB1 ity=0\n\
             \x20 max rel err    : 0.000e0 at grid#0\n\
             \x20 ulp histogram  : 0:7950\n",
        )];
        let points = error_points(&runs);
        assert_eq!(points.len(), 2);
        let svg = strip(1, 1e-13, &points).expect("a figure");
        let text = svg.finish("test fixture");
        // One data dot plus the legend swatch.
        assert_eq!(text.matches(r#"class="bad""#).count(), 2, "{text}");
        assert!(text.contains("exactly 0"), "{text}");
        assert_eq!(
            outside_the_gate(&points, 1e-13),
            vec!["tier1 TESUB1 ity=2, offset as f64".to_string()]
        );
        assert!(outside_the_gate(&points, 1.0).is_empty());
    }

    /// The positive control contributes 7,950 comparisons at `>=16` ULP, which
    /// would otherwise be the only thing the histogram showed.
    #[test]
    fn the_ulp_figure_leaves_out_the_blocks_beyond_the_gate() {
        let runs = [run(
            "tier1_enthalpy",
            Libm::Vendored,
            "test whatever ... tier1 TESUB1 ity=2, offset as f64\n\
             \x20 max rel err    : 1.597e-5 at grid#704\n\
             \x20 ulp histogram  : >=16:7950\n\
             tier1 TESUB1 ity=0\n\
             \x20 max rel err    : 0.000e0 at grid#0\n\
             \x20 ulp histogram  : 0:7950\n",
        )];
        let svg = ulp(1, 1e-13, &runs).expect("a figure");
        let text = svg.finish("test fixture");
        assert!(
            text.contains("7,950 of 7,950 identical to the last bit"),
            "{text}"
        );
        assert!(!text.contains("&gt;=16 ULP"), "{text}");
        assert!(
            text.contains("One block beyond the gate is left out"),
            "{text}"
        );
    }

    #[test]
    fn the_tier_four_sweep_is_read_row_by_row() {
        let runs = [
            run(
                "tier4_trajectory",
                Libm::Vendored,
                "test tier4_fault_scenarios ... Tier 4, fault scenarios, 4 h each:\n\
                 \x20 nominal    within XNS 4.000 h, first out none, worst 2.45e-5\n\
                 \x20 IDV(1)     within XNS 4.000 h, first out none, worst 1.60e-9\n",
            ),
            run(
                "tier4_trajectory",
                Libm::Platform,
                "\x20 nominal    within XNS 4.000 h, first out none, worst 0.00e0\n",
            ),
        ];
        let points = noise_points(&runs);
        assert_eq!(points.len(), 3);
        assert_eq!(points[0].scenario, "nominal");
        assert!((points[0].hours - 4.0).abs() < 1e-12);
        assert!((points[1].worst - 1.60e-9).abs() < 1e-20);
        assert_eq!(points[2].libm, Libm::Platform);
        assert_eq!(points[2].worst, 0.0);

        let svg = noise_band(&points, 4.0).expect("a figure");
        let text = svg.finish("test fixture");
        assert!(text.contains("1 of 3 runs are bit-identical"), "{text}");
        assert!(text.contains("class=\"altr\""), "{text}");
    }

    /// The nominal 8 h test prints its hours and its error on separate lines,
    /// which is a different horizon and must not join the 4 h sweep.
    #[test]
    fn the_eight_hour_block_is_not_mistaken_for_a_sweep_row() {
        let runs = [run(
            "tier4_trajectory",
            Libm::Vendored,
            "test tier4_nominal_trajectory ... Tier 4, nominal, vendored libm, 8 h:\n\
             \x20 within XNS for : 8.000 h\n\
             \x20 worst at end   : 5.149e-5 of XNS\n",
        )];
        assert!(noise_points(&runs).is_empty());
    }

    #[test]
    fn the_calibration_reads_both_line_shapes() {
        let runs = [run(
            "tier5_battery",
            Libm::Vendored,
            "test a_stuck_valve ... XMV(10) under IDV(14): sd 6.9, shifted 69.4\n\
             \x20 KS cross 1.0000 vs within Some(0.075)\n\
             === nominal (4 seeds, reporting only) ===\n\
             \x20 frobenius: cross 6.973021e-12, within-source max 1.090611e1, p 1.0000\n",
        )];
        let points = calibration_points(&runs);
        assert_eq!(points.len(), 2);
        assert_eq!(points[0].what, "KS, one variable shifted 10 sd");
        assert!((points[0].cross - 1.0).abs() < 1e-12);
        assert_eq!(points[1].what, "correlation matrix, nominal");
        assert!(points[1].cross < points[1].within);

        let svg = calibration(&points, "smoke battery").expect("a figure");
        let text = svg.finish("test fixture");
        // Each point is coloured by its own two values and by nothing else:
        // the shifted variable is on the wrong side of the diagonal, the
        // correlation matrix is on the right one.
        assert!(
            text.contains(
                r#"class="bad" cx="625.0" cy="179.1" r="5.5"><title>KS, one variable shifted"#
            ),
            "{text}"
        );
        assert!(
            text.contains(r#"class="ok" cx="346.1" cy="125.1" r="5.5"><title>correlation matrix"#),
            "{text}"
        );
        assert!(text.contains("1 of 2 points lie above it"), "{text}");
        assert!(!text.contains("unlabelled for want of room"), "{text}");
        // Both axes must have the same pixels per decade, or the diagonal is
        // not the line of equality and the whole figure is a geometric lie.
        assert!(
            text.contains(r#"<line class="gate" x1="300.0" y1="476.0" x2="700.0" y2="76.0"/>"#)
        );
    }

    /// A caption lives inside a raw HTML block, where mdBook parses no
    /// markdown at all, so its markup has to be converted rather than passed
    /// through looking like backticks.
    /// mdBook copies these files to the site as well as inlining them, and a
    /// figure opened on its own is parsed as XML, where a comment containing
    /// `--` is fatal. Every command that writes one has `--tiers` in it, so
    /// this is not a hypothetical. A processing instruction avoids that and
    /// makes mdBook's own HTML parser log a parse error instead, which is how
    /// the provenance ended up in `<metadata>`.
    #[test]
    fn a_figures_provenance_survives_both_parsers() {
        let root = std::env::temp_dir().join(format!("xtask-plot-{}", std::process::id()));
        let _ = std::fs::remove_dir_all(&root);
        let runs = [run(
            "tier2_heat",
            Libm::Vendored,
            "test t ... tier1 heat QUS\n\x20 max rel err : 9.572e-13 at x\n",
        )];
        let svg = strip(2, 1e-12, &error_points(&runs)).expect("a figure");
        write_figure(&root, &svg, "cargo xtask validate --tiers 1,2,3 --smoke").expect("write");

        let path = root.join(format!("{DIR}/tier2-errors.svg"));
        let text = std::fs::read_to_string(&path).expect("read back");
        let mut lines = text.lines();
        assert!(
            lines.next().expect("a first line").starts_with("<svg id="),
            "{text}"
        );
        let second = lines.next().expect("a second line");
        assert!(second.starts_with("<metadata>GENERATED by "), "{second}");
        assert!(second.contains("--tiers 1,2,3 --smoke"), "{second}");
        assert!(!text.contains("<!--"), "an XML comment survived: {text}");
        assert!(
            !text.contains("<?"),
            "a processing instruction survived: {text}"
        );

        // And the index can still read it back off the file.
        let (by, commit) =
            crate::report::read_provenance(&root, &format!("{DIR}/tier2-errors.svg"))
                .expect("provenance");
        assert_eq!(by, "cargo xtask validate --tiers 1,2,3 --smoke");
        assert!(!commit.is_empty());
        let _ = std::fs::remove_dir_all(&root);
    }

    #[test]
    fn a_caption_becomes_inline_html() {
        assert_eq!(
            inline_html("the `1e-12` gate & the **band**"),
            "the <code>1e-12</code> gate &amp; the <strong>band</strong>"
        );
        assert_eq!(inline_html("a < b"), "a &lt; b");
        // An unpaired marker is text, not the start of a bold run that eats
        // the rest of the caption.
        assert_eq!(inline_html("2 ** 3"), "2 ** 3");
    }

    /// Tier 3's situation: everything bit-identical, nothing to put on a log
    /// axis. The page says so instead of drawing an empty one.
    #[test]
    fn a_tier_that_differed_nowhere_gets_no_error_axis() {
        let runs = [run(
            "tier1_disturbance",
            Libm::Vendored,
            "test t ... tier1 TESUB5 ADIST\n\
             \x20 max rel err : 0.000e0 at seed#0\n\
             \x20 ulp histogram : 0:20000\n",
        )];
        let points = error_points(&runs);
        assert_eq!(points.len(), 1);
        assert!(strip(3, 1e-13, &points).is_none());
        // The claim still gets made, in the units that can carry it.
        assert!(ulp(3, 1e-13, &runs).is_some());
    }

    /// The full battery is twenty-one scenarios, and their correlation-matrix
    /// points land within a pixel of one another. Stacking twenty-one labels
    /// would bury the figure, so the ones with no room are dropped to their
    /// tooltips and the figure says how many.
    #[test]
    fn a_crowded_calibration_drops_the_labels_it_cannot_place() {
        let points: Vec<CalibrationPoint> = (0..21)
            .map(|i| CalibrationPoint {
                what: format!("correlation matrix, IDV({i})"),
                cross: 4e-12,
                within: 1.1e1,
            })
            .collect();
        let svg = calibration(&points, "full battery").expect("a figure");
        let text = svg.finish("test fixture");
        assert!(text.contains("unlabelled for want of room"), "{text}");
        // Every point still has its dot and its tooltip, whatever happened to
        // its label.
        assert_eq!(
            text.matches("<title>correlation matrix, IDV(").count(),
            21,
            "{text}"
        );
    }

    #[test]
    fn an_empty_run_draws_no_figure() {
        assert!(strip(2, 1e-12, &[]).is_none());
        assert!(noise_band(&[], 4.0).is_none());
        assert!(calibration(&[], "none").is_none());
        assert!(ulp(2, 1e-12, &[]).is_none());
    }

    #[test]
    fn a_count_is_grouped_for_reading() {
        assert_eq!(group(9_987_490), "9,987,490");
        assert_eq!(group(0), "0");
        assert_eq!(group(999), "999");
        assert_eq!(group(1000), "1,000");
    }

    /// Every figure must survive being pasted into a page: no raw `<` or `&`
    /// from a label, and one `id` so two figures on a page cannot collide.
    #[test]
    fn labels_are_escaped_into_the_document() {
        let runs = [run(
            "tier2_balances",
            Libm::Vendored,
            "test t ... YP(1..50) & the <gate>\n\x20 max rel err : 6.093e-14 at x\n",
        )];
        let svg = strip(2, 1e-12, &error_points(&runs)).expect("a figure");
        let text = svg.finish("test fixture");
        assert!(text.contains("&amp; the &lt;gate&gt;"), "{text}");
        assert!(!text.contains("& the <gate>"), "{text}");
        assert!(text.contains(r#"id="tier2-errors""#), "{text}");
    }
}
