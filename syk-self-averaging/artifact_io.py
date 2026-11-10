"""Shared artifact writing for the campaign producers.

write_json_atomic dumps to a sibling temp file and then os.replace()s it
over the target, so a kill mid-dump can never leave a truncated artifact
behind.  The campaigns that use it run from minutes to most of a day.
The bytes written are exactly those of json.dump(payload, fh,
indent=indent), so switching a producer to it does not change its
artifact format.
"""

from __future__ import annotations

import json
import os
from pathlib import Path


def write_json_atomic(path, payload, indent=1) -> None:
    path = Path(path)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(payload, fh, indent=indent)
    os.replace(tmp, path)
