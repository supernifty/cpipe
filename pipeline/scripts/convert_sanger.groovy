convert_sanger = {
        filter("phred") {
            exec """
                    gunzip -c  $input.gz | maq sol2sanger  - - | gzip -c > $output.gz
            """
        }
}

run {
   "%.gz" * [ convert_sanger ]
}
