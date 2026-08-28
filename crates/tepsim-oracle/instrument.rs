//! Build-time instrumentation of the original Fortran.
//!
//! Included by `build.rs`. Not part of the crate's own source tree.
//!
//! # Why this exists
//!
//! The original keeps `ISD`, the shutdown flag, as a local variable in
//! `TEFUNC`, so there is no symbol to link against and a differential test
//! cannot see whether the plant tripped. Upstream's f2py variant solves this by
//! editing the source; we cannot, because `reference/` is ground truth and is
//! checksummed on every gate run.
//!
//! So the edit happens at build time, into `OUT_DIR`, and the vendored file is
//! never touched.
//!
//! # Why string replacement rather than a patch file
//!
//! Applying a `.patch` needs the `patch` binary and fails in confusing ways
//! when context drifts. Explicit replacements that *assert their own
//! pre-image* fail with a precise message naming the text that was expected and
//! not found, which is what a future session actually needs. Each edit below
//! also has to match exactly once: matching twice would mean the anchor is not
//! as specific as it looks.

/// One textual edit, with the exact text it expects to find.
struct Edit {
    /// What this edit is for, quoted back in the panic message.
    what: &'static str,
    /// Must occur exactly once in the source.
    from: &'static str,
    to: &'static str,
}

/// Hoist `ISD` out of `TEFUNC`'s locals and into `COMMON/SHUTDN/`.
///
/// This mirrors what upstream's f2py variant does, and is behaviour-preserving
/// for the same reason: `TEFUNC` assigns `ISD=0` unconditionally at
/// `teprob.f:702`, before the eight shutdown tests and before both reads, so
/// the separate borrow of `ISD` as a `TESUB7` sign flag at `teprob.f:387`
/// cannot leak into the shutdown decision. Splitting that borrow into its own
/// `IRAND` local removes the overlap entirely.
///
/// A test asserts that the instrumented build agrees with the pristine one on
/// the derivative vector, so "behaviour-preserving" is checked, not asserted.
const EDITS: &[Edit] = &[
    Edit {
        what: "declare ISD in COMMON/SHUTDN/ within TEFUNC, and add an IRAND local",
        from: "      INTEGER NN,I,ISD\n",
        to: "      INTEGER ISD\n      COMMON/SHUTDN/ISD\n      INTEGER NN,I,IRAND\n",
    },
    Edit {
        what: "use IRAND, not ISD, as the TESUB7 sign flag in the DO 910 loop",
        from: "      ISD=-1\n      HWLK=HSPAN(I)*TESUB7(ISD)+HZERO(I)\n",
        to: "      IRAND=-1\n      HWLK=HSPAN(I)*TESUB7(IRAND)+HZERO(I)\n",
    },
    Edit {
        what: "declare ISD in COMMON/SHUTDN/ within TEINIT",
        from: "      DOUBLE PRECISION XMEAS,XMV\n      COMMON/PV/XMEAS(41),XMV(12)\n      INTEGER IDV\n      COMMON/DVEC/IDV(20)\n      DOUBLE PRECISION G\n      COMMON/RANDSD/G\n",
        to: "      DOUBLE PRECISION XMEAS,XMV\n      COMMON/PV/XMEAS(41),XMV(12)\n      INTEGER IDV\n      COMMON/DVEC/IDV(20)\n      DOUBLE PRECISION G\n      COMMON/RANDSD/G\n      INTEGER ISD\n      COMMON/SHUTDN/ISD\n",
    },
    // Anchored on the DO 550 terminator, not on `TNEXT(I)=0.1D0`: that line
    // also appears in TEFUNC's t=0 reset block, so it matches twice. The
    // uniqueness assertion below caught that before it could do any damage.
    Edit {
        what: "initialise ISD to zero at the end of TEINIT",
        from: "  550 CONTINUE\n      TIME=0.0\n",
        to: "  550 CONTINUE\n      ISD=0\n      TIME=0.0\n",
    },
    // Tier 3. Every draw in the file comes through `TESUB7`, so one edit here
    // records the whole stream. It is appended *after* both assignments to the
    // function result, so nothing it does can affect the value returned; the
    // pristine-probe test in `tests/instrumentation_equivalence.rs` checks
    // that rather than leaving it as an argument.
    //
    // `TRCN` counts every draw, but only the first `TRCCAP` are stored. So an
    // overflow is visible to the reader as a count past the capacity rather
    // than as a silently truncated trace.
    //
    // It *saturates* at `TRCCAP+1` rather than counting on. `TRCN` is a
    // Fortran `INTEGER`, which is 32 bits, and nothing clears it between the
    // evaluations of a long run: Tier 5's battery makes about 264 draws per
    // step over 172,800 steps per run, so it passes 2^31 after roughly fifty
    // runs. The counter then goes *negative*, `TRCN.LE.TRCCAP` becomes true
    // again, and `TRCVAL(TRCN)` writes outside the array. That is a segfault
    // deep inside the Fortran, tens of millions of steps into a battery, with
    // nothing in the immediate vicinity to suggest why. B-0047b found it that
    // way. Saturating costs one comparison and makes the failure impossible.
    //
    // The sign flag is recorded as well as the value, because `TESUB7` returns
    // two different scalings depending on it (`teprob.f:1552-1553`), and a
    // port could call the wrong one the right number of times.
    Edit {
        what: "record every TESUB7 draw into COMMON/RNGTRC/ for Tier 3",
        from: "      COMMON/RANDSD/G\n      G=DMOD(G*9228907.D0,4294967296.D0)\n      IF(I.GE.0)TESUB7=G/4294967296.D0\n      IF(I.LT.0)TESUB7=2.D0*G/4294967296.D0-1.D0\n",
        to: "      COMMON/RANDSD/G\n      INTEGER TRCN,TRCSGN,TRCCAP\n      PARAMETER (TRCCAP=4096)\n      DOUBLE PRECISION TRCVAL\n      COMMON/RNGTRC/TRCVAL(TRCCAP),TRCSGN(TRCCAP),TRCN\n      G=DMOD(G*9228907.D0,4294967296.D0)\n      IF(I.GE.0)TESUB7=G/4294967296.D0\n      IF(I.LT.0)TESUB7=2.D0*G/4294967296.D0-1.D0\n      TRCN=TRCN+1\n      IF(TRCN.LE.TRCCAP)THEN\n      TRCVAL(TRCN)=TESUB7\n      TRCSGN(TRCN)=1\n      IF(I.LT.0)TRCSGN(TRCN)=-1\n      ELSE\n      TRCN=TRCCAP+1\n      ENDIF\n",
    },
];

/// The one edit `temain_mod.f` needs.
///
/// The file is an *unnamed main program*: it has no `PROGRAM` statement, so
/// gfortran compiles it to `MAIN__` and emits a `main` that calls it. Linking
/// that into a Rust test binary would collide with the harness's own entry
/// point, and worse, running it would open fifteen files under `~/` and
/// simulate 48 hours.
///
/// Turning the main program into a subroutine nothing calls solves both. The
/// nineteen `CONTRLn` subroutines below it are unaffected and become linkable,
/// which is the whole point: they are what Phase 4 ports.
///
/// The anchor is the one comment in the file with *two* spaces before
/// "MEASUREMENT"; the twenty-three copies inside the subroutines have three.
/// That is exactly the kind of near-miss the uniqueness assertion exists for.
const DRIVER_EDITS: &[Edit] = &[Edit {
    what: "make temain_mod.f's unnamed main program a subroutine nothing calls",
    from: "C  MEASUREMENT AND VALVE COMMON BLOCK\nC\n      DOUBLE PRECISION XMEAS, XMV\n",
    to: "C  MEASUREMENT AND VALVE COMMON BLOCK\nC\n      SUBROUTINE TEMAIN_UNUSED\n      DOUBLE PRECISION XMEAS, XMV\n",
}];

/// Apply the driver edit. See [`DRIVER_EDITS`].
pub(crate) fn instrument_driver(source: &str) -> String {
    apply(source, DRIVER_EDITS)
}

/// Apply every edit, or panic with a message that says which one failed and why.
pub(crate) fn instrument(source: &str) -> String {
    apply(source, EDITS)
}

fn apply(source: &str, edits: &[Edit]) -> String {
    let mut out = source.to_string();
    for edit in edits {
        let hits = out.matches(edit.from).count();
        assert!(
            hits != 0,
            "instrumentation failed: could not find the text to {what}.\n\
             The vendored Fortran no longer contains:\n{from:?}\n\
             This means reference/fortran/teprob.f changed. Do not loosen this \
             match; work out what changed and why.",
            what = edit.what,
            from = edit.from,
        );
        assert!(
            hits == 1,
            "instrumentation failed: the anchor to {what} matched {hits} times, \
             expected exactly once. The anchor is not specific enough, and \
             applying it would edit the wrong places.",
            what = edit.what,
        );
        out = out.replace(edit.from, edit.to);
    }
    out
}
