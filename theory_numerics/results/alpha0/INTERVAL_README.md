# Interval verification of exact protected-response overlaps

Run from the repository root after installing \`requirements.txt\`:

\`\`\`bash
python theory_numerics/verify_overlap_sign_intervals.py
\`\`\`

The script reads the 33 published \`explicit_overlaps_alpha0.csv\` cases and
computes directed interval enclosures for the unnormalized **positive-half-axis**
endpoint charge integral at \`t=1\`. It splits \`[0,14 sigma]\` into 1000 intervals
at 35 decimal-digit interval precision. A polynomial Gaussian-tail bound encloses
the omitted infinite range. For every case, the test compares the interval upper
endpoint with zero **before** converting any number to a float.

The closest upper enclosure is about -2.6882447e-6 at \`(sigma,L)=(0.15,8)\`.
All 33 upper enclosures are negative. The checked-in CSV records rounded diagnostic
upper enclosures and tail bounds. Running the script generates the full interval
endpoint strings in \`theory_numerics/reproduced_alpha0/interval_endpoint_signs_alpha0.json\`.

This validates signs for the specified one-dimensional alpha=0 protected response
and combines with the intermediate value theorem to give protected-coordinate
overlaps. It does not prove equivalence of full QFT outcome distributions.
