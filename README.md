After knitting and rendering `bird_summary.Rmd`, overwrite this file with a symlink to `bird_summary.md`.

E.g., on Linux/macOS from terminal, in the directory that contains the module:

```bash
cd bird_summary
rm README.md && ln -s bird_summary.md README.md
```
