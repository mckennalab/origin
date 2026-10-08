# README figures

Both images are rendered by the code blocks in the top-level `README.md`,
verbatim and with `seed = 1`, so they show what that code actually produces
rather than a tidied version of it.

To regenerate after changing an example, extract the blocks and run them:

```bash
Rscript -e '
  s <- readLines("README.md")
  starts <- grep("^```r$", s); ends <- grep("^```$", s)
  blocks <- Map(function(a, b) s[(a + 1):(b - 1)],
                starts, ends[sapply(starts, function(a) which(ends > a)[1])])
  blocks <- Filter(function(b) !any(grepl("install.packages", b)), blocks)
  code <- unlist(blocks); code <- code[!startsWith(code, "#>")]
  i <- max(grep("^plot\\\\(tree", code))
  png("docs/images/quickstart_lineage.png", 1000, 1000, res = 130)
  eval(parse(text = paste(code[1:i], collapse = "\n"))); dev.off()
  png("docs/images/edits_on_tree.png", 1500, 1100, res = 130)
  eval(parse(text = paste(code[(i + 1):length(code)], collapse = "\n"))); dev.off()
'
```
