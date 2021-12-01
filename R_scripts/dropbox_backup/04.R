# npp
plj <- scan("File4.txt")

######
my.file <- file("File4.txt")
lines <- readLines(my.file)
close(my.file)
######

cat("xyz\nsafa\n", file = "test_file.txt")
### encoding, UTF8, UTF16, \n, \t, \r
cat("XXX")
print("XXX")

######
fl1 <- read.table("File1.txt", header = T, sep = "\t");
fl2 <- read.csv("File2.txt", header = T, sep = "\t");


######

reg.exmpl <- c("krava", "plava", "prava", "tava", "vatra", "bla", "blabla", "blablabla", 
               "abcabcabc", "abcblaabc", "blaXbla", "abcYabc", "abc1z", "abc4z")

reg.exmpl

grep("va", reg.exmpl)

grep("va", reg.exmpl, value=T)

grep("^va", reg.exmpl)

grep("ra$", reg.exmpl, value=T)

grep("[lt]ava", reg.exmpl, value=T)

grep("bla", reg.exmpl, value=T)

grep("(bla){3}", reg.exmpl, value=T)

grep("(abc|bla){3}", reg.exmpl, value=T)

grep("(bla)?", reg.exmpl, value=T)

grep("(bla).", reg.exmpl, value=T)

grep("(bla)*", reg.exmpl, value=T)

grep("p.ava", reg.exmpl, value=T)

grep("p[^r]ava", reg.exmpl, value=T)

grep("(abc|bla)[XY]\\1", reg.exmpl, value=T)

grep("abc\\dz", reg.exmpl, value=T)

# regexone.com