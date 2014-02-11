convert_sanger = {
        var MAQ : "maq"

        filter("phred") {
            exec """
                    gunzip -c  $input.gz | $MAQ sol2sanger  - - | gzip -c > $output.gz
            """
        }
}

run {
   "%.gz" * [ convert_sanger ]
}
