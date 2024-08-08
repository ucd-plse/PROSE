import sys, os

sys.path.insert(
    0,
    os.path.realpath(
        os.path.join(
            os.path.join(
                os.path.dirname(
                    os.path.realpath(
                        __file__
                    )
                ),
                os.pardir
            ),
            os.path.join(
                "src",
                "python"
            )
        )
    )
)

from transformers import ProseProjectTransformer
from precimonious import PrecimoniousSearch
from bruteforce import BruteForceSearch
from profiling import preprocess_project

del sys.path[0], sys, os