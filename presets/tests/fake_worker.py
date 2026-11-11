"""Stand-in for recipe/_align_worker.py (no kalign): copies input to output.
A case whose input path contains "hang" sleeps forever, "crash" exits with 7, "cpu" burns ~0.3 s CPU."""
import json
import os
import shutil
import sys
import time


def do(inp, out):
    c0 = time.process_time()
    if "hang" in inp:
        time.sleep(3600)
    if "crash" in inp:
        os._exit(7)
    if "cpu" in inp:
        t = time.process_time()
        while time.process_time() - t < 0.3:
            pass
    shutil.copy(inp, out)
    return {"wall": 0.01, "cpu": time.process_time() - c0, "error": None}


if sys.argv[1] == "--serve":
    proto = os.fdopen(os.dup(1), "w", buffering=1)
    os.dup2(2, 1)
    print("noise on stdout", flush=True)  # must not corrupt the protocol
    proto.write(json.dumps({"ready": True}) + "\n")
    for line in sys.stdin:
        req = json.loads(line)
        proto.write(json.dumps(do(req["in"], req["out"])) + "\n")
else:
    job = json.load(open(sys.argv[1]))
    rows = [do(i, o) for i, o in job["cases"]]
    json.dump({"rows": rows}, open(sys.argv[2], "w"))
