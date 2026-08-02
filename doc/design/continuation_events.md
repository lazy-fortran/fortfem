---
title: Continuation event diagnostics
---

# Continuation event diagnostics

`classify_continuation_event` compares caller-owned signed margins at two
successive continuation points. It reports the first deterministic sign
crossing, or the first margin whose absolute value is within a caller-supplied
tolerance:

\[
  m_i^{old}m_i^{new}<0 \Longrightarrow \text{sign crossing},
  \qquad
  \min(|m_i^{old}|,|m_i^{new}|)\le\varepsilon
  \Longrightarrow \text{near zero}.
\]

The minimum margin over all entries is returned for step-size control. The
same neutral diagnostic can represent a level-set cut, separatrix, resonance,
interface, or branch-selection margin; the caller supplies the physical
meaning. Event classification is discrete by design. A client must freeze or
invalidate JVP/VJP products when a topology event is reported rather than
silently differentiating across a changed algebraic space.

## API

```fortran
call classify_continuation_event( &
    previous_margin, current_margin, tolerance, event_code, event_index, &
    minimum_margin, status)
```

The public codes are `CONTINUATION_EVENT_NONE`,
`CONTINUATION_EVENT_SIGN_CROSSING`, and `CONTINUATION_EVENT_NEAR_ZERO`.
Inputs must have equal, non-empty shapes and finite values; the tolerance must
be non-negative and finite. The independent test covers crossing, near-zero,
no-event, deterministic index selection, and invalid tolerance cases.
