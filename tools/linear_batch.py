"""Python transport for TIAGO's normalized linear Gaussian batch executable."""
from pathlib import Path
import subprocess
import tempfile

import numpy as np


def infer(executable, cases):
    """Each case supplies whitened design/data, normalization and prior SDs."""
    with tempfile.TemporaryDirectory() as directory:
        source, target = Path(directory)/'input.txt', Path(directory)/'output.txt'
        with source.open('w') as stream:
            stream.write(f'{len(cases)}\n')
            for case in cases:
                design, data = np.asarray(case['design']), np.asarray(case['data'])
                m, n = design.shape
                stream.write(f"{m} {n} {case['log_normalization']:.17g}\n")
                np.savetxt(stream, np.asarray(case['prior_sd'])[None], fmt='%.17g')
                np.savetxt(stream, np.c_[data, design], fmt='%.17g')
        subprocess.run([str(executable), str(source), str(target)], check=True)
        results = []
        with target.open() as stream:
            for index in range(len(cases)):
                number, count, evidence, loss = stream.readline().split()
                if int(number) != index+1:
                    raise ValueError('unexpected linear posterior index')
                count = int(count)
                mean = np.fromstring(stream.readline(), sep=' ')
                covariance = np.array([np.fromstring(stream.readline(), sep=' ')
                                       for _ in range(count)])
                if mean.shape != (count,) or covariance.shape != (count, count):
                    raise ValueError('incomplete linear posterior')
                results.append(dict(mean=mean, covariance=covariance,
                                    log_evidence=float(evidence), loss=float(loss)))
    return results
