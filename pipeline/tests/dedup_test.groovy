
load '../config.groovy'

// All the core pipeline stages in the pipeline
load '../pipeline_stages_config.groovy'

run {
  dedup
}
