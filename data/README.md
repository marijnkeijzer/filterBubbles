# Analysis data

Each row is one participant-topic observation. `study1.csv` contains 312 rows (156 per topic); `study2.csv` contains 474 rows (237 per topic). Respondents contributed to two topics, but respondent identifiers and linkage have been removed. Rows were shuffled independently within each released file. Original values and missingness are preserved for retained columns; no noise, imputation, or synthetic observations were added. Missing values are `NA`.

| Variable | Meaning / coding |
| --- | --- |
| `topic` | `dev`: development aid; `tax`: taxation topic. |
| `t0_gen` | Pre-stimulus general opinion, rescaled to [0,1]. |
| `t1_gen` | Immediate post-stimulus general opinion, [0,1]. |
| `t2_gen` | Follow-up general opinion, [0,1]; Study 1 only. |
| `sti_auth` | Respondent's perceived opinion of the stimulus author, [0,1]. |
| `leftright` | Respondent's left-right ideological identification, 0–100. |
| `sti_lr` | Respondent's perceived left-right identification of the stimulus author, stored on a 0–1 scale. The original model's transformation is preserved in code/common.R. |
| `ExpCondition` | Original numeric assignment code; interpret jointly with study/topic. Study 2 analysis marks exposure as `ExpCondition != 0` for dev and `ExpCondition != 1` for tax. |
| `stimulus_direction` | Original stimulus direction, 0 contra / 1 pro. |
| `stimulus_morality` | Original binary moral-frame code (0/1); Study 1 only. Retained as coded in the original model. |

The exported Study 1 data have 10 columns and Study 2 has 8. Model-specific exclusions are applied in `code/common.R`:

- Main Study 1 models: complete `t0_gen`, `t1_gen`, `sti_auth` (311 rows).
- Reference and lagged Study 1 models: additionally complete `t2_gen` (124 rows).
- Study 2 models: complete `t0_gen`, `t1_gen`, `sti_auth`, `sti_lr`, `leftright`, then the topic-specific exposure filter above (306 rows).
- Descriptive table: complete `t0_gen`, `t1_gen`, `sti_auth`, `sti_lr`, `leftright`, following the original descriptive script. This is distinct from the model sample filter.

Starting from model-ready data means consent checks, survey recoding and cleaning decisions precede this release. Do not reinterpret the exported number of rows as the recruited participant count. The pooled models treat topic rows as observations; respondent-clustered alternatives cannot be fitted without a linkage key, which is intentionally absent.

## Posterior draws

`posterior/` contains the original parameter draws needed for model predictions and ABM calibration. They contain no participant records: only model parameters, chain numbers and retained-draw indices. These are replication inputs, not validation reports.
