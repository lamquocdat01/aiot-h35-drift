"""H35: check rendered PDF - pages, abstract words, cited DOIs printed, '??', placeholders.
Usage: python check_render.py <main.pdf> <main.bbl>"""
import re, sys, json
import pymupdf

doc = pymupdf.open(sys.argv[1]); text = "".join(p.get_text() for p in doc)
bbl = open(sys.argv[2], encoding="utf-8").read()
items = bbl.split(chr(92) + "bibitem")[1:]
dois = sorted(set(d.lower().rstrip(".") for d in re.findall(r"(10\.\d{4,9}/[^\s}{]+)", bbl)))
flat = re.sub(r"\s+", "", text).lower().replace("‐", "-")
missing = [d for d in dois if d not in flat]
ab = text.split("Abstract", 1)[1].split("Keywords", 1)[0]
kw = text.split("Keywords:", 1)[1].split("\n1", 1)[0] if "Keywords:" in text else ""
res = {"pages": doc.page_count, "abstract_words": len(re.findall(r"[A-Za-z0-9][\w\-.%']*", ab)),
       "keywords": len([k for k in kw.split(";") if k.strip()]), "cited": len(items), "dois_in_bbl": len(dois),
       "dois_missing_in_pdf": missing, "question_marks": text.count("??"), "placeholders": re.findall(r"<[A-Z_]+>", text)}
print(json.dumps(res, indent=1))
