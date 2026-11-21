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
open(sys.argv[2], "w", encoding="utf-8").write(s)
print("bib ok", s.count("howpublished"))
