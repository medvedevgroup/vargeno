# Algorithm Description
## Author: Chen Sun(chensunx@gmail.com)
The encoded k-mer is as follows:
1. For SNP k-mer and Lite Ref k-mer
----$hi_{24}$----|--------$lo_{40}$--------
2. For Ref k-mer
------$hi_{32}$------|------$lo_{32}$------
We use $hi$ and $lo$ for general case.

$BF$ stores $lo_{40}$ information, if query $lo_{40}$ does not exist , then do not consider any high changes.

For general case, a k-mer can be divided into three parts, illustrated as follows:
----$k_{41-63}$----|---$k_{32-40}$---|--------$k_{0-31}$--------

k-mers can be divided into three categories:
1. SNP k-mers, SK
2. Ref k-mers, RK
3. Lite Ref k-mers, L-RK

For $k_{0-31}$, for all k-mers, if block-size < threshold, then linear direct match; if block-size >= threshold, then binary search.

For $k_{32-40}$, 
1. For SK, if block-size < threshold, then linear direct match; if block-size >= threshold, then binary search.
2. For RK, all go to binary search.
3. For L-RK, if block-size < threshold, then linear direct match; if block-size >= threshold, then binary search.

For $k_{41-63}$, for all k-mers, binary search.

However, if we write like this, then divide the code into 2+5+1=8 separated code blocks. 

Here to save the code blocks in VarGeno, we implementated as follows:

### Code Block 1
```
if block-size >= threshold,
	binary search k0 to k31
else
	call linear direct match function about lo part
```
Here about $lo$ part
1. For SK, it's $lo_{40}$
2. For RK, it's $lo_{32}$
3. For L-RK, it's $lo_{40}$

Then we need to deal with kmer hi part.

Note here also need to check if $lo$ part is in BF.  If not, we do not even need to iterate $hi$ part.

Since there are three types of k-mer, and most of them need binary search.
1. For SK, if block-size >= threshold, from k32 to k63, binary search. Else, from k40 to k63, binary search.
2. For RK, if block-size >= threshold, from k32 to k63, binary search. Else, from k32 to k63, binary search. 
3. For L-RK, if block-size >= threshold, from k32 to k63, binary search. Else, from k40 to k63, binary search.

### Code Block 2
```
for i in 32 to 63,
	// for ref k-mer
	if i >= 40, or block-size >= threshold, or ref is lite:
		if lo part in BF:
			bianry search ref k-mer k_i
	
	// for snp k-mer
	if i >= 40 or block-size >= threshold:
		if lo part in BF:
			binary search snp k-mer k_i
```
Here also need to check if $lo$ part is in BF.  If not, we do not even need to iterate $hi$ part.