This is an attempt to improve the CDK HashFingerprint (Fingerprinter class).
The idea behind the improved version is borrowed from my blog improvised hashing function and their impact on the fingerprints. 

http://chembioinfo.com/2011/10/30/revisiting-molecular-hashed-fingerprints/

**Compile with dependencies**

````
 mvn clean compile assembly:single
````

## Command line interface ##

**Test improved CDK FP**

````
 java -jar target/fingerprinter-1.0-SNAPSHOT-jar-with-dependencies.jar mol scaffold 250
````   

**Improved CDK HashedFingerprinter class with 1024 size FP**


CASES |	TP	| FP |	TN	| FN	|	ACCURACY	|	TPR	|	FPR	|	Time (mins)
------|----|----|----|----|----------|-----|-----|-------------
25*25	|	25	|1	|597|	2|	0.996|		0.926|	0.002|	0
50*50	|	51	|9|	2433	|7	|0.994|		0.880	|0.004|	0
75*75		|76	|21|	5502|	26	|0.992|		0.746	|0.004	|0
100*100	|	101	|47	|9784|	68|	0.989	|	0.598	|0.005	|0.01
125*125		|129	|71|	15337	|88|	0.990|		0.595	|0.005|	0.01
150*150		|154	|79|	22155	|112|	0.992	|	0.579	|0.004|	0.01
175*175	|	183	|106	|30070	|266	|0.988		|0.408	|0.004	|0.01
200*200	|	210	|137	|39330	|323|	0.989		|0.394|	0.004|	0.02
225*225	|	236	|149	|49875|	365	|0.990		|0.393	|0.003|	0.02
250*250		|266|	225	|61489	|520	|0.989|		0.339|	0.004|	0.02

**Test CDK default FP**
 
````
 java -jar target/fingerprinter-1.0-SNAPSHOT-jar-with-dependencies.jar mol cdk 250
````

**CDK Default Fingerprinter class with 1024 size FP**


CASES |	TP	| FP |	TN	| FN	|	ACCURACY	|	TPR	|	FPR	|	Time (mins)
------|----|----|----|----|----------|-----|-----|-------------
25*25	|	25	|4	|594|	2|	0.991|		0.926|	0.007|	0
50*50	|	51	|20|	2422|	7|	0.990	|	0.880 |	0.009|	0
75*75		|76|	68|	5455|	26|	0.984|		0.746	|0.013|	0
100*100	|	101|	181	|9650|	68|	0.976		|0.598	|0.019|	0.01
125*125	|	129	|257	|15151	|88	|0.978		|0.595|	0.017	|0.01
150*150	|	154|	325	|21909	|112	|0.981		|0.579|	0.015|	0.01
175*175	|	183|	648	|29528	|266	|0.971		|0.408|	0.022|	0.01
200*200	|	210|	810	|38657	|323	|0.972		|0.394	|0.021	|0.02
225*225	|	236|	928	|49096	|365	|0.975		|0.393	|0.019	|0.02
250*250	|	266	|1240	|60474	|520|	0.972|		0.339	|0.021|	0.03


**Note: New Scaffold fingerprinter reduces number of False Positives (FP) - High Accuracy**
