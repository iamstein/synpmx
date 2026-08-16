## This is the author's To do list, to be maintained by him

- Try with a couple more datasets (ask at next dose finding toolbox meeting)

- Spot check the code (ask co-author and colleague if they'd like to help here)
    - consider /ultrareview
    - validating scorecard (spot checks)
    - validating avatar algorithm (spot checks)
    
- In parallel - Meet with an internal colleague (mid august)

    Four that decide something.  Do not leave without these:
    - where is synadam output allowed to go?  off GxP, laptop, cloud tenant,
      external LLM API?  Decides whether the precedent covers our use case.
    - what did you have to produce?  The artifact list is the spec for the
      spot-check work, so it decides how deep that goes and when.
    - can synpmx ride along on the precedent, or is it its own path from zero?
      Decides the timeline.
    - does a new release trigger re-approval?  If yes, freeze synpmx sooner.

    Context:
    - what exactly got approved - the package as software, a process, or each
      dataset release?  blanket or per-study sign-off?  conditions or expiry?
    - how long did it take, who signed off, one person or a committee?
    - what got pushed back on, and what surprised you?
    - was there a formal risk assessment?
    - what level of software validation was required?
    - did anyone independently read the code, and did that have to be documented?
    - did you have to demonstrate anything adversarially?
    - did the reviewers object at the outset to whole trajectories?
    - what features were most important (e.g. .rds?  anything else)
    - who owns and maintains it now?
    - would you co-sponsor, or review the privacy argument before I take it up?

- Make poster (due Sept 28 for online submission)
    - Find the external publication clearance path and its lead time.  Separate
      from tool approval, and the long pole if it is slow.  Ask the colleague.
    - The abstract says the synthetic data does not contain real patient data.
      AVATAR blends from real patient values, so the poster has to correct that
      rather than repeat it.  Abstract was written before the method existed.

## Broader GenAI activities
- Move WRITING_FOR_ANDY.md to TrinityMetrics/skills/ so it can serve
  more than one project, with AGENTS.md pointing synpmx at its raw URL instead of
  keeping a copy. On the way out, qualify the commit citations as synpmx@c972373,
  and leave the contract table and standing caveat here as synpmx-specific.