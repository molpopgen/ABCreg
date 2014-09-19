##ABCreg

The citation for this software is doi:10.1186/1471-2156-10-35, and can be found [here](http://www.biomedcentral.com/1471-2156/10/35).

Release 0.1.0 corresponds to the publised version.

The current (master) version differs from the published in that the output files are not gzipped using [zlib](http://zlib.net).

An addition to the documentation is that users of bash-like shells may use gzipped input files via implicit subshells:

```{sh}
reg -p <(gunzip -c prior.gz) -d <(gunzip -c data.gz) -P 1 -S 1 -b output -T -t 0.001
```
