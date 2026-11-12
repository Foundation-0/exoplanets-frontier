# Plain-Language Guide: What This Project Does and What We Found

This guide explains the work in `Foundation-0/exoplanets-frontier` in everyday terms, without specialist astronomy jargon, so anyone on the team or in the public can see what has been built, what failed, what succeeded, and why it matters.

## 1. The Big Picture

NASA's TESS space telescope stares at millions of stars watching for tiny, repeating dips in starlight caused when an orbiting planet crosses in front of its star (a "transit"). Earlier work in [eumemic/exoplanets](https://github.com/eumemic/exoplanets) scanned over 4 million stars in archival TESS images and flagged **665 candidate planets** that were not yet in official NASA catalogs.

Turning a raw brightness dip into a real, follow-up-ready world requires solving three hard practical problems:

1. **Clock drift over years**: If a planet's orbit period is off by even 5 seconds in 2018, after 1,000 orbits by 2026 the planet arrives more than an hour early or late. Looking only at the old predicted minute makes the planet seem to have vanished, and aiming a major observatory like the James Webb Space Telescope at the wrong hour wastes precious telescope time.
2. **Blurry camera pixels ("photobombers")**: Each pixel on TESS's camera covers a large patch of sky (21 arcseconds across, dozens of times blurrier than ground telescopes). Often, a faint pair of eclipsing stars sitting next door inside the same blurry pixel dims deeply and mimics a small planet on the main star.
3. **Starspots on small red stars**: Small, cool red stars ("M-dwarfs") are the best places to find rocky, Earth-sized planets, because a small planet blocks a bigger fraction of a small star. However, many red stars spin fast and are covered in dark starspots that make their brightness rise and fall like a big ocean wave, hiding tiny planetary dips underneath.

## 2. What Stalled at First, and How Pulling the Threads Unraveled Breakthroughs

In our first campaign, we tested held-out 2026 TESS data (Sectors 106 and 107) and hunted for new planets around cool red stars:
- **22 candidate planets showed up right on schedule** in fresh 2026 observations (`results/DISCOVERIES.md`).
- Several other candidates looked too faint or shifted against the edge of our 60-minute clock window, two candidates previously dismissed as "likely nearby false positives" unexpectedly reappeared right on schedule, and **all 19 of our red-star triggers turned out to be spinning starspots** rather than planets.

Instead of abandoning those stalled runs, we built targeted fixes (`despot.py`, `ttv_track.py`, `blend_gaia.py`, `catalog_665.py`, `converge.py`) and ran them on live NASA and European Space Agency (Gaia) star catalogs:

1. **Rescuing "lost" planets by widening the clock window (2018 to 2026)**:
   Expanding our timing tracker from 1 hour to 4 hours across 8 years of TESS observations proved that several downgraded candidates had simply drifted slightly ahead or behind their 2018 clocks. Fixing a 44-second period offset on candidate **`TIC 166848769`** more than doubled its signal strength (from 6.3 to 13.5 times the noise level), and fixing a 5.2-second offset rescued **`TIC 332644900`** from apparent non-detection while revealing signs of possible gravitational tugs from another body. Meanwhile, an apparent timing wobble we had flagged earlier on **`TIC 70557792`** resolved cleanly into a simple 1.45-second clock refinement.
2. **Pinpointing which star inside a blurry pixel actually dimmed**:
   By subtracting in-transit images from out-of-transit images and mapping every neighboring star from the Gaia space astrometry catalog (accounting for how stars move across the sky between 2016 and 2026), we localized the dimming on **`TIC 253923533`** (previously written off as a nearby false alarm) to within a quarter of a pixel (`4.9 arcseconds`) of the target star itself. We also built a catalog-wide filter across all **665 `eumemic/exoplanets` candidates** (`results/eumemic_gaia_dr4_readiness.csv`) ahead of the European Space Agency's 25 November 2026 archive shutdown, identifying "Platinum" targets where **zero** neighboring stars in the pixel are bright enough to fake the planet dip.
3. **Peeling away starspot waves without erasing planets**:
   We built a two-pass mathematical wave remover (`harmonic_despot`) that fits and subtracts a star's spinning spot pattern while ignoring sharp downward planet dips. Tested blindly on a known spotted red star (**TOI-540**), naive software got fooled by the 17.4-hour starspot wave, whereas our wave remover stripped the spots and recovered the sub-Earth rocky planet **TOI-540 b** to within 1.9 seconds of its published orbit, with zero false alarms on scrambled data.

## 3. Where the Highest-Value Outputs Live

- **`results/platinum_followup_targets.csv` & `results/eumemic_gaia_dr4_readiness.csv`**: Ranked, dilution-checked catalog of the 665 `eumemic/exoplanets` candidates, highlighting clean solo stars, temperate (liquid-water-zone) worlds, and top targets for atmospheric study with the James Webb Space Telescope, pre-indexed by Gaia ID for the 2 December 2026 Gaia DR4 release.
- **`results/ttv_ephemeris_audit.csv` & `results/blend_localization.csv`**: Updated 2018-to-2026 orbital clocks and sub-pixel localization measurements so observatories know which candidates are on-target and the exact minute their future transits occur.
- **`results/despotted_mdwarf_audit.csv` & `results/vetted_new_candidates.csv`**: Full ledger of our two-pass starspot-cleaned cool-dwarf planet search, recording both any passing candidates and every vetted negative result so future runs never repeat finished stars.
