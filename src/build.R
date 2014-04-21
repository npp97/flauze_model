require(compiler)
system('R CMD SHLIB flauz.c');
cmpfile('flauz_subs.r');

