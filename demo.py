"""Three short examples on the cubic path."""
from monotone_history import decode


def examples():
    zero = [["0", "0", "0"], ["0", "0", "0"]]
    return {
        "unique": {
            "times": ["-.02", "0"],
            "positions": [["-.0396", "-.019404", "0"], ["0", "0", "0"]],
            "linf_errors": [".0005", ".0005"],
        },
        "ambiguous": {
            "times": ["-.01", ".01"],
            "positions": zero,
            "linf_errors": [".06", ".06"],
        },
        "incompatible": {
            "times": ["-.01", ".01"],
            "positions": zero,
            "linf_errors": [".0005", ".0005"],
        },
    }


if __name__ == "__main__":
    for name, history in examples().items():
        result = decode(history, "1/2")
        if result["status"] != name.upper():
            raise RuntimeError(f"Unexpected result for {name}: {result['status']}")
        branch = result["selected"] or "none"
        print(f"{name}: {result['status']}; selected branch: {branch}")
