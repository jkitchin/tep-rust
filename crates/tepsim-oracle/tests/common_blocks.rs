//! Verifies the Rust mirrors of the Fortran `COMMON` blocks are laid out
//! correctly, and that the build-time instrumentation changes no numbers.
//!
//! A wrong offset here would not crash. It would silently return a neighbouring
//! variable, and every differential test downstream would compare the wrong
//! quantities while looking perfectly healthy. So the layout is checked by
//! reading fields whose values are known from `TEINIT`, spread across the whole
//! block including the far end, rather than by trusting the struct definition.

#![cfg(feature = "oracle")]

use std::mem::size_of;

use tepsim_oracle::{Const, Oracle, Teproc, Wlk};

/// Exact float comparison, by bits, so the intent is unambiguous and clippy's
/// `float_cmp` has nothing to object to. These are stored constants, not
/// computed results: they either round-trip exactly or the layout is wrong.
#[track_caller]
fn assert_exact(actual: f64, expected: f64, what: &str) {
    assert_eq!(
        actual.to_bits(),
        expected.to_bits(),
        "{what}: got {actual:?}, expected {expected:?}"
    );
}

/// The value a Fortran literal *without* a `D` suffix actually ends up with.
///
/// In fixed-form Fortran, `25.4` is a single-precision constant. Assigning it
/// to a `DOUBLE PRECISION` variable widens the already-rounded `f32`, so the
/// stored value is not the nearest double to 25.4. 182 of the constant
/// assignments in `TEINIT` are written this way; see the test at the bottom of
/// this file, and B-0006.
fn single(literal: f32) -> f64 {
    f64::from(literal)
}

#[test]
fn block_sizes_match_the_fortran_declarations() {
    assert_eq!(
        size_of::<Teproc>(),
        Teproc::LEN_BYTES,
        "/TEPROC/ is 580 doubles then 12 integers; a mismatch means a field was \
         miscounted or the compiler inserted padding"
    );
    assert_eq!(size_of::<Wlk>(), Wlk::LEN_BYTES);
    assert_eq!(size_of::<Const>(), Const::LEN_BYTES);
    assert_eq!(Teproc::LEN_BYTES, 4688);
    assert_eq!(Wlk::LEN_BYTES, 1104);
    assert_eq!(Const::LEN_BYTES, 896);
}

/// Molecular weights, `teprob.f:941-948`. If `/CONST/` were misaligned these
/// would come back as Antoine coefficients instead, which look nothing alike.
#[test]
fn const_block_holds_the_molecular_weights() {
    let mut oracle = Oracle::lock();
    oracle.init();
    let c = oracle.constants();
    // XMW is the last of fourteen arrays in /CONST/, so reading it correctly
    // means every earlier array is the right length too.
    //
    // Note 25.4: `teprob.f:942` writes it without a D suffix, so the stored
    // value is the f32 rounding, not the nearest double. Every other weight
    // here happens to be exactly representable in both.
    let expected = [2.0, single(25.4), 28.0, 32.0, 46.0, 48.0, 62.0, 76.0];
    for (i, (got, want)) in c.xmw.iter().zip(expected).enumerate() {
        assert_exact(*got, want, &format!("XMW({})", i + 1));
    }
    assert_ne!(
        c.xmw[1].to_bits(),
        25.4_f64.to_bits(),
        "XMW(2) must NOT be the nearest double to 25.4. If this ever passes, \
         the Fortran was changed, and a Rust port using 25.4_f64 would be \
         wrong by 1.5e-8 relative: five orders of magnitude past the Tier 1 \
         tolerance."
    );
}

/// Vessel volumes and coolant holdups, early and middle in `/TEPROC/`.
/// From `teprob.f:1118-1125`.
#[test]
fn teproc_holds_the_vessel_constants() {
    let mut oracle = Oracle::lock();
    oracle.init();
    let p = oracle.teproc();

    // These are exactly representable in f32 and f64 alike, so the missing D
    // suffix costs nothing.
    assert_exact(p.vtr, 1300.0, "VTR reactor total volume");
    assert_exact(p.vts, 3500.0, "VTS separator total volume");
    assert_exact(p.vtc, 156.5, "VTC stripper total volume");
    assert_exact(p.vtv, 5000.0, "VTV mixing zone total volume");
    assert_exact(p.hwr, 7060.0, "HWR reactor coolant holdup");
    assert_exact(p.hws, 11138.0, "HWS condenser coolant holdup");
    assert_exact(p.cpflmx, 280275.0, "CPFLMX compressor max flow");

    // Written with a D suffix at teprob.f:1122-1123, so full double precision.
    assert_exact(p.htr[0], 0.06899381054, "HTR(1) heat of reaction 1");
    assert_exact(p.htr[1], 0.05, "HTR(2) heat of reaction 2");

    // Written without one at teprob.f:1171, so f32 precision.
    assert_exact(
        p.cpprmx,
        single(1.3),
        "CPPRMX compressor max pressure ratio",
    );
}

/// `XNS` sits at offset 525 of 580, near the far end of `/TEPROC/`. Reading it
/// correctly means the whole run of preceding fields is right.
/// From `teprob.f:1256-1296`.
#[test]
fn teproc_holds_the_measurement_noise_scales() {
    let mut oracle = Oracle::lock();
    oracle.init();
    let p = oracle.teproc();

    // These carry D suffixes at teprob.f:1256-1296, so they are exact doubles.
    assert_exact(p.xns[0], 0.0012, "XNS(1), A feed flow noise");
    assert_exact(p.xns[1], 18.0, "XNS(2), D feed flow noise");
    assert_exact(p.xns[40], 0.5, "XNS(41), the very last noise scale");
    assert!(
        p.xns.iter().all(|v| *v > 0.0),
        "every measurement has a positive noise scale"
    );
}

/// `VST` and `IVST` are the final two fields, and `IVST` is the only integer
/// array in the block. If the double-to-integer transition were misplaced,
/// `IVST` would read as garbage reinterpreted from a double's bit pattern.
#[test]
fn teproc_ends_with_vst_then_the_ivst_integers() {
    let mut oracle = Oracle::lock();
    oracle.init();
    let p = oracle.teproc();

    for (i, v) in p.vst.iter().enumerate() {
        assert_exact(*v, 2.0, &format!("VST({})", i + 1));
    }
    assert_eq!(
        p.ivst, [0; 12],
        "IVST(I)=0 for all twelve, teprob.f:1107. A misplaced boundary here \
         would show up as huge or nonsensical integers."
    );
}

/// `/WLK/` ends the same way: eleven double arrays then twelve integers.
/// `TEINIT` sets every `TNEXT` to 0.1, which is also why nothing is drawn at t=0.
#[test]
fn wlk_block_is_initialised_as_teinit_leaves_it() {
    let mut oracle = Oracle::lock();
    oracle.init();
    let w = oracle.wlk();

    for (i, v) in w.tnext.iter().enumerate() {
        assert_exact(*v, 0.1, &format!("TNEXT({}), teprob.f:1362", i + 1));
    }
    for (i, v) in w.tlast.iter().enumerate() {
        assert_exact(*v, 0.0, &format!("TLAST({})", i + 1));
    }
    assert_exact(w.hspan[11], 1.5, "HSPAN(12), teprob.f:1355");
    assert_exact(w.hzero[11], 2.0, "HZERO(12)");
    assert_eq!(
        w.idvwlk, [0; 12],
        "the trailing integer array, zero with no disturbances active"
    );
}

/// Fortran stores `XST(8,13)` column-major, so `XST(i,j)` is `xst[j-1][i-1]`.
/// Transposing this would still typecheck and still be the right size, so it
/// is checked against a composition that has to sum to one.
#[test]
fn two_dimensional_arrays_are_indexed_stream_then_component() {
    let mut oracle = Oracle::lock();
    let (t, yy) = oracle.init();
    oracle.derivatives(t, &yy);
    let p = oracle.teproc();

    // Stream 4 is the A/C feed, whose three components are set every call at
    // teprob.f:407-410 and must sum to one.
    // Tolerance is 1e-7, not 1e-12, and that is not slack: the feed
    // compositions at teprob.f:1134-1159 are written as literals like 0.9999
    // and 0.0001 with no D suffix, so they are f32 values widened to double and
    // sum to 1 only to about 1.7e-8. A tighter bound here would be asserting
    // something the original does not do.
    for stream in [0usize, 1, 2, 3] {
        let total: f64 = p.xst[stream].iter().sum();
        assert!(
            (total - 1.0).abs() < 1e-7,
            "stream {} composition sums to {total}, not 1. If this is far off \
             rather than slightly off, the array is being read transposed.",
            stream + 1
        );
    }
}

/// The shutdown flag exists only because `build.rs` hoists `ISD` into a
/// `COMMON` block. At the nominal operating point the plant has not tripped.
#[test]
fn the_shutdown_flag_is_reachable_and_clear_at_nominal() {
    let mut oracle = Oracle::lock();
    let (t, yy) = oracle.init();
    oracle.derivatives(t, &yy);
    assert_eq!(
        oracle.shutdown_flag(),
        0,
        "the nominal operating point must not be a shutdown condition"
    );
}

/// Drive the plant into a shutdown condition and confirm the flag reports it.
///
/// `teprob.f:706` trips on reactor temperature above 175 C. Rather than
/// contrive a state, this raises reactor temperature directly through the
/// energy holdup.
#[test]
fn the_shutdown_flag_trips_on_reactor_overtemperature() {
    let mut oracle = Oracle::lock();
    let (t, yy) = oracle.init();

    // Sanity: not tripped to begin with.
    oracle.derivatives(t, &yy);
    assert_eq!(oracle.shutdown_flag(), 0);

    // YY(9) is reactor internal energy. Raising it lifts reactor temperature.
    let mut hot = yy;
    hot[8] *= 1.5;
    let yp = oracle.derivatives(1.0 / 3600.0, &hot);

    assert_eq!(
        oracle.shutdown_flag(),
        1,
        "reactor temperature should have exceeded the 175 C trip at teprob.f:706"
    );
    assert!(
        yp.iter().all(|v| *v == 0.0),
        "on a trip the original zeroes all 50 derivatives at teprob.f:807, \
         freezing the plant rather than shutting it down. That is the Class C \
         quirk the delta register has to decide about."
    );
}

/// The whole `/TEPROC/` block round-trips: write a modified snapshot, read it
/// back byte-identical.
#[test]
fn teproc_round_trips_through_the_common_block() {
    let mut oracle = Oracle::lock();
    oracle.init();

    let original = oracle.teproc();
    let mut modified = original;
    modified.tcr = 123.456;
    modified.xns[40] = 0.75;
    modified.ivst[11] = 1;
    modified.fcm[12][7] = -2.5;

    oracle.set_teproc(&modified);
    let read_back = oracle.teproc();
    assert_eq!(read_back, modified, "/TEPROC/ must round-trip exactly");
    assert_ne!(
        read_back, original,
        "and the write must actually have landed"
    );

    oracle.set_teproc(&original);
    assert_eq!(oracle.teproc(), original);
}

#[test]
fn wlk_round_trips_through_the_common_block() {
    let mut oracle = Oracle::lock();
    oracle.init();

    let original = oracle.wlk();
    let mut modified = original;
    modified.adist[0] = 0.25;
    modified.tnext[11] = 9.5;
    modified.idvwlk[6] = 1;

    oracle.set_wlk(&modified);
    assert_eq!(oracle.wlk(), modified);

    oracle.set_wlk(&original);
    assert_eq!(oracle.wlk(), original);
}

/// The single most dangerous transcription hazard in the whole port.
///
/// Fixed-form Fortran treats a literal with no exponent letter, or with an `E`
/// exponent, as *single* precision. Assigning it to a `DOUBLE PRECISION`
/// variable widens a value that has already been rounded to 24 bits of
/// mantissa. So `XMW(2)=25.4` at `teprob.f:942` does not store the nearest
/// double to 25.4; it stores 25.399999618530273.
///
/// The error is about 1.5e-8 relative. Tier 1 asks for 1e-13. A Rust port that
/// transcribes these constants as ordinary `f64` literals therefore fails by
/// five orders of magnitude, and would look like a deep numerical problem
/// rather than a transcription one.
///
/// B-0006 must decide per constant, from the presence or absence of a `D`
/// suffix in the source, which of the two values to use.
#[test]
fn constants_without_a_d_suffix_carry_only_single_precision() {
    let mut oracle = Oracle::lock();
    oracle.init();
    let c = oracle.constants();

    let stored = c.xmw[1];
    let naive = 25.4_f64;
    let faithful = f64::from(25.4_f32);

    assert_exact(stored, faithful, "XMW(2) as the Fortran actually stores it");
    assert_ne!(stored.to_bits(), naive.to_bits());

    let relative = (stored - naive).abs() / naive;
    assert!(
        relative > 1e-9 && relative < 1e-7,
        "the f32 widening error on XMW(2) should be about 1.5e-8, got {relative:e}"
    );
}
