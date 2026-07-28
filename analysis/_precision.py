"""Serialisation precision for saved statistics.

Byte-level drift gates are the cheapest way to catch a pipeline that has
silently changed. They are defeated by floating-point noise: SciPy versions
differ in the last few digits of a p-value (1.7868222549680278e-23 vs
...208e-23) without any change to the science. Serialising at a documented
precision makes the gate portable across environments while still catching
any difference that could matter.

12 significant figures is far beyond what is reported or interpretable, and
well inside the digits that vary between library versions.
"""

SIGFIGS = 12


def round_sig(x, sig=SIGFIGS):
    """Round a float to `sig` significant figures. Passes through non-floats."""
    import math
    if not isinstance(x, float) or x == 0 or not math.isfinite(x):
        return x
    return round(x, -int(math.floor(math.log10(abs(x)))) + (sig - 1))


def clean(obj, sig=SIGFIGS):
    """Recursively round every float in a nested structure before writing."""
    if isinstance(obj, dict):
        return {k: clean(v, sig) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [clean(v, sig) for v in obj]
    return round_sig(obj, sig)
