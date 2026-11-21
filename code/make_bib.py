"""H35: lit/refs.bib -> manuscript/refs.bib for unsrtnat: @misc entries get howpublished with the DOI URL,
typographic hyphen/quote normalised to ASCII. Usage: python make_bib.py <in.bib> <out.bib>"""
import re, sys
s = open(sys.argv[1], encoding="utf-8").read()
s = s.replace("‐", "-").replace("‑", "-").replace("’", "'").replace("‘", "'")


def fix(m):
    body = m.group(0)
    d = re.search(r"doi\s*=\s*\{([^}]+)\}", body)
    if body.lower().startswith("@misc") and d:
        url = "\\url{https://doi.org/%s}" % d.group(1)
        h = re.search(r"howpublished\s*=\s*\{([^}]*)\}", body)
        if h:
            body = body.replace(h.group(0), "howpublished = {%s, %s}" % (h.group(1), url))
        else:
            body = body.rstrip()[:-1].rstrip().rstrip(",") + ",\n  howpublished = {%s}\n}" % url
    return body


s = re.sub(r"@\w+\{[^@]*?\n\}", fix, s, flags=re.S)

# F3/H5: protect capitals from the style's lower-casing, fix HTML entities and DataCite author split; drop "none" issue numbers.
FIXES = [
    ("M\\&amp;V cascade for LNG", "{M\\&V} cascade for {LNG}"),
    ("in Korea's national", "in {Korea}'s national"),
    ("{W-UDDR: A Unified Framework for Automated Drift Detection in IoT-based Air Quality Monitoring Systems}",
     "{{W-UDDR}: A Unified Framework for Automated Drift Detection in {IoT}-based Air Quality Monitoring Systems}"),
    ("Time-uniform Chernoff bounds", "Time-uniform {Chernoff} bounds"),
    ("delay time in page's sequential", "delay time in {Page}'s sequential"),
    ("{Page's sequential procedure", "{{Page}'s sequential procedure"),
    ("{SEQUENTIAL MONITORING OF CHANGES IN DYNAMIC LINEAR MODELS, APPLIED TO THE U.S. HOUSING MARKET}",
     "{Sequential monitoring of changes in dynamic linear models, applied to the {U.S.} housing market}"),
    ("author = {Vito, Saverio}", "author = {De Vito, Saverio}"),
]
for a, b in FIXES:
    assert a in s, a
    s = s.replace(a, b)
s = re.sub(r"\n\s*number = \{none\},", "", s)
open(sys.argv[2], "w", encoding="utf-8").write(s)
print("bib ok", s.count("howpublished"))
