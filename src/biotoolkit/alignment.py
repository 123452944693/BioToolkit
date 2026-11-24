"""Needleman-Wunsch global DNA alignment with a linear gap penalty."""
from dataclasses import dataclass
from .sequence import normalize_dna

@dataclass(frozen=True)
class Alignment:
    first: str
    second: str
    score: int

def global_align(first: str, second: str, *, match: int = 2,
                 mismatch: int = -1, gap: int = -2) -> Alignment:
    """Return one optimal alignment; tie order diagonal, up, left.

    IUPAC symbols are compared literally. Quadratic memory/time; limit
    matrices to one million cells. Scores must be integers, gap negative.
    """
    if not all(isinstance(s, int) for s in (match, mismatch, gap)) or gap >= 0:
        raise ValueError("Integer scores and a negative gap penalty required")
    a, b = normalize_dna(first), normalize_dna(second)
    m, n = len(a), len(b)
    if (m+1)*(n+1) > 1_000_000:
        raise ValueError("Alignment exceeds the one-million-cell limit")
    dp = [[0]*(n+1) for _ in range(m+1)]
    for i in range(1,m+1): dp[i][0] = i*gap
    for j in range(1,n+1): dp[0][j] = j*gap
    for i in range(1,m+1):
        for j in range(1,n+1):
            diag = dp[i-1][j-1] + (match if a[i-1] == b[j-1] else mismatch)
            dp[i][j] = max(diag, dp[i-1][j]+gap, dp[i][j-1]+gap)
    x, y = [], []
    i, j = m, n
    while i or j:
        if i and j and dp[i][j] == dp[i-1][j-1] + (match if a[i-1] == b[j-1] else mismatch):
            x.append(a[i-1]); y.append(b[j-1]); i -= 1; j -= 1
        elif i and dp[i][j] == dp[i-1][j] + gap:
            x.append(a[i-1]); y.append("-"); i -= 1
        else:
            x.append("-"); y.append(b[j-1]); j -= 1
    return Alignment("".join(reversed(x)), "".join(reversed(y)), dp[m][n])
