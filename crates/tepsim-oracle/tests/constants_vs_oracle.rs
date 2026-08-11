//! Compares the transcribed constant table against the values gfortran actually
//! stored in `COMMON/CONST/`.
//!
//! This is the decisive check. `tepsim-core`'s own test re-derives the values by
//! parsing the source, which shares an assumption with the transcription: that
//! my reading of Fortran literal semantics is right. This one asks the compiler
//! instead, so it would catch a misunderstanding that both of the others share.

#![cfg(feature = "oracle")]

use tepsim_core::constants;
use tepsim_oracle::Oracle;

#[test]
fn every_constant_matches_what_gfortran_stored() {
    let mut oracle = Oracle::lock();
    oracle.init();
    let c = oracle.constants();

    let pairs: [(&str, &[f64; 8], &[f64; 8]); 14] = [
        ("AVP", constants::AVP.as_array(), &c.avp),
        ("BVP", constants::BVP.as_array(), &c.bvp),
        ("CVP", constants::CVP.as_array(), &c.cvp),
        ("AH", constants::AH.as_array(), &c.ah),
        ("BH", constants::BH.as_array(), &c.bh),
        ("CH", constants::CH.as_array(), &c.ch),
        ("AG", constants::AG.as_array(), &c.ag),
        ("BG", constants::BG.as_array(), &c.bg),
        ("CG", constants::CG.as_array(), &c.cg),
        ("AV", constants::AV.as_array(), &c.av),
        ("AD", constants::AD.as_array(), &c.ad),
        ("BD", constants::BD.as_array(), &c.bd),
        ("CD", constants::CD.as_array(), &c.cd),
        ("XMW", constants::XMW.as_array(), &c.xmw),
    ];

    let mut wrong = Vec::new();
    let mut checked = 0;
    for (name, ours, theirs) in pairs {
        for i in 0..8 {
            checked += 1;
            if ours[i].to_bits() != theirs[i].to_bits() {
                wrong.push(format!(
                    "  {}({}): Rust {:?} vs Fortran {:?}  (differ by {:.3e} relative)",
                    name,
                    i + 1,
                    ours[i],
                    theirs[i],
                    if theirs[i] == 0.0 {
                        (ours[i] - theirs[i]).abs()
                    } else {
                        (ours[i] - theirs[i]).abs() / theirs[i].abs()
                    }
                ));
            }
        }
    }

    assert_eq!(checked, 112, "all 112 constants must be compared");
    assert!(
        wrong.is_empty(),
        "{} of 112 constants differ from what gfortran stored:\n{}\n\nThe usual \
         cause is transcribing a literal without its D suffix as a plain f64, or \
         wrapping a D-suffixed one in single().",
        wrong.len(),
        wrong.join("\n")
    );
}

/// Guards against the comparison passing because both sides are zero-filled.
#[test]
fn the_oracle_block_is_actually_populated() {
    let mut oracle = Oracle::lock();
    oracle.init();
    let c = oracle.constants();
    let nonzero = [&c.avp, &c.ah, &c.ag, &c.ad, &c.xmw]
        .iter()
        .flat_map(|a| a.iter())
        .filter(|v| **v != 0.0)
        .count();
    assert!(
        nonzero > 20,
        "only {nonzero} non-zero values in a sample of COMMON/CONST/; TEINIT \
         cannot have run, so a comparison against it would prove nothing"
    );
}
