"""Export small CSV tables for the C++ example."""
import json
import os

ROUND = r"D:\CodingPrograms\Protein Design\round"
FIN = os.path.join(ROUND, "result", "fin")
PEAKDIR = os.path.join(ROUND, "result", "peak", "fin_peaks")

with open(os.path.join(FIN, "dictionary.json"), encoding="utf-8") as f:
    dictionary = json.load(f)

with open(os.path.join(FIN, "dictionary.csv"), "w", encoding="utf-8") as f:
    f.write("center,id,phi,psi\n")
    for center, peaks in dictionary.items():
        for p in peaks:
            f.write(f"{center},{p['id']},{p['phi']},{p['psi']}\n")

with open(os.path.join(FIN, "extras.csv"), "w", encoding="utf-8") as f:
    f.write("a,b,c,phi,psi,weight\n")
    with open(os.path.join(FIN, "contexts.jsonl"), encoding="utf-8") as g:
        for line in g:
            rec = json.loads(line)
            for e in rec["extra"]:
                f.write(f"{rec['a']},{rec['b']},{rec['c']},{e['phi']},{e['psi']},{e['weight']}\n")

print("csv written")
