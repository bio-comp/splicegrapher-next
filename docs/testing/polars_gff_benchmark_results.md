# Polars GFF Single-Cycle Evaluation

## Decision
- recommendation: `defer`
- decision: `threshold_not_met`
- evaluated real datasets: `3`
- compared real datasets: `3`
- winning datasets: `0`
- thresholds: runtime>=20%, memory<=+10%
- notes: required runtime and memory win thresholds not met

## Synthetic
### small
#### ingest
| dataset | loader | status | rows | mean_s | max_s | peak_mib |
|---|---|---:|---:|---:|---:|---:|
| small | gene_model | ok | 500 | 0.019650 | 0.025381 | 0.553 |
| small | rows | ok | 500 | 0.003431 | 0.004310 | 0.407 |
| small | polars_df | ok | 500 | 0.003964 | 0.004979 | 0.407 |

#### end_to_end
| dataset | loader | status | rows | mean_s | max_s | peak_mib |
|---|---|---:|---:|---:|---:|---:|
| small | gene_model | ok | 300 | 0.018961 | 0.019426 | 0.645 |
| small | rows | ok | 300 | 0.007098 | 0.008441 | 0.407 |
| small | polars_df | ok | 300 | 0.008121 | 0.008185 | 0.411 |

### medium
#### ingest
| dataset | loader | status | rows | mean_s | max_s | peak_mib |
|---|---|---:|---:|---:|---:|---:|
| medium | gene_model | ok | 3750 | 0.136482 | 0.141784 | 4.258 |
| medium | rows | ok | 3750 | 0.022392 | 0.022468 | 3.006 |
| medium | polars_df | ok | 3750 | 0.026144 | 0.027002 | 3.006 |

#### end_to_end
| dataset | loader | status | rows | mean_s | max_s | peak_mib |
|---|---|---:|---:|---:|---:|---:|
| medium | gene_model | ok | 2250 | 0.156574 | 0.170964 | 5.297 |
| medium | rows | ok | 2250 | 0.048310 | 0.049439 | 3.029 |
| medium | polars_df | ok | 2250 | 0.063136 | 0.063695 | 3.294 |

### large
#### ingest
| dataset | loader | status | rows | mean_s | max_s | peak_mib |
|---|---|---:|---:|---:|---:|---:|
| large | gene_model | ok | 15000 | 0.557787 | 0.590906 | 17.282 |
| large | rows | ok | 15000 | 0.090897 | 0.091280 | 12.035 |
| large | polars_df | ok | 15000 | 0.108949 | 0.110545 | 12.035 |

#### end_to_end
| dataset | loader | status | rows | mean_s | max_s | peak_mib |
|---|---|---:|---:|---:|---:|---:|
| large | gene_model | ok | 9000 | 0.660991 | 0.679746 | 22.779 |
| large | rows | ok | 9000 | 0.217922 | 0.233263 | 12.627 |
| large | polars_df | ok | 9000 | 0.277665 | 0.283683 | 14.206 |

## Real
### at5
#### ingest
| dataset | loader | status | rows | mean_s | max_s | peak_mib |
|---|---|---:|---:|---:|---:|---:|
| at5 | gene_model | ok | 27 | 0.000956 | 0.000974 | 0.036 |
| at5 | rows | ok | 27 | 0.000174 | 0.000178 | 0.033 |
| at5 | polars_df | ok | 27 | 0.000221 | 0.000228 | 0.033 |

#### end_to_end
| dataset | loader | status | rows | mean_s | max_s | peak_mib |
|---|---|---:|---:|---:|---:|---:|
| at5 | gene_model | ok | 11 | 0.001090 | 0.001095 | 0.036 |
| at5 | rows | ok | 11 | 0.000317 | 0.000331 | 0.033 |
| at5 | polars_df | ok | 11 | 0.000480 | 0.000489 | 0.033 |

### at1
#### ingest
| dataset | loader | status | rows | mean_s | max_s | peak_mib |
|---|---|---:|---:|---:|---:|---:|
| at1 | gene_model | ok | 53 | 0.001748 | 0.001800 | 0.042 |
| at1 | rows | ok | 53 | 0.000301 | 0.000307 | 0.051 |
| at1 | polars_df | ok | 53 | 0.000381 | 0.000399 | 0.051 |

#### end_to_end
| dataset | loader | status | rows | mean_s | max_s | peak_mib |
|---|---|---:|---:|---:|---:|---:|
| at1 | gene_model | ok | 22 | 0.001977 | 0.001989 | 0.042 |
| at1 | rows | ok | 22 | 0.000611 | 0.000637 | 0.051 |
| at1 | polars_df | ok | 22 | 0.000770 | 0.000777 | 0.051 |

### at2
#### ingest
| dataset | loader | status | rows | mean_s | max_s | peak_mib |
|---|---|---:|---:|---:|---:|---:|
| at2 | gene_model | error:ValueError | 8 | 0.000000 | 0.000000 | 0.000 |
| at2 | rows | ok | 8 | 0.000092 | 0.000097 | 0.020 |
| at2 | polars_df | ok | 8 | 0.000125 | 0.000128 | 0.020 |

#### end_to_end
| dataset | loader | status | rows | mean_s | max_s | peak_mib |
|---|---|---:|---:|---:|---:|---:|
| at2 | gene_model | error:ValueError | 8 | 0.000000 | 0.000000 | 0.000 |
| at2 | rows | ok | 0 | 0.000106 | 0.000110 | 0.020 |
| at2 | polars_df | ok | 0 | 0.000171 | 0.000175 | 0.020 |
