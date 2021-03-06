# NIP (Network Interaction Parser)
In order to avoid the arguments curse, we split the NIP into two parts, the data preprocessing and parser. 

##  Installation
Download NIP package and install from the directory (Python version : 3.x):
```
git clone https://github.com/bruce803/NIP.git
cd  NIP
```

# Get it started:
## Reading Data

Before the network analysis, we read all raw data into some temporal python pickle files [https://docs.python.org/3/library/pickle.html].  Generally, we have four kinds of data, two gene lists, one global network, and pathway/function data. This process can be done with the “readData” command. 
Suppose that global network(net.txt), in gene list (cle.txt), out gene list (ple.txt), pathway (pathway_list, a directory) are stored in the current directory, one example of the dataRead command works as follows,
Example：
```
python readData.py -n net-karate.txt -i in.txt -o out.txt -f  pathway
```
The readData step returns us four python pickle files storing in the directory “.\data”. Actually, we do not have to care about these temporal files. It is just the preparation work for the next step—network parser.


## Network parser

In this section, we introduce the details of implementation of NIP toolbox.  The input variables includes a public gene network, two gene lists. Our target is to search the interactions between two gene lists over the public network. Four operations (intersect, merge, filter, delete) are provided for the convenience of achieving that goal.
The overview of this tool box is,

![A demo of NIP](https://github.com/bruce803/NIP/blob/master/result/demo-NIP.png)



### interact
As illustrated by the Figure below, the public network was divided into three parts. The green and cyan lists are the two significant gene sets that we interested in, the red cluster (a-e) is the bridge between them.  The four green (x1, x2, x3, x4) and cyan nodes (y1, y2, y3, y4) are the shared genes between the public network and the two gene lists. In fact, the total number of the green or cyan nodes is more than 4. Figure 1 only plots the intersected part with the global network. Our target is to find out how the green list interacts with cyan list via the public network.  In this paper, we only focus on three different cases. The first case is that the green and cyan nodes may directly connected. That is the one hop case. If the number of directed connection is zero, we can further search the two hops links (just as the picture shows). Our toolbox can at most support the three hops case, namely the green and cyan nodes are bridged by two layers nodes.
These three cases corresponding to the three choices of the “-l” and ”--layer” arguments. For example, the following command running for the third case (three hops),

```Python run.py -l 2```

The default setting is layer=2.

1 hop             |  2 hops
:-------------------------:|:-------------------------:
![A 1 hop example.](https://github.com/bruce803/NIP/blob/master/result/intersection0.png)  |  ![A 2 hops example.](https://github.com/bruce803/NIP/blob/master/result/interaction2hops.png)


###	merge
Merge the nodes in the middle layers. Usually, there are too many nodes in the intersected layer, so we may prefer to merge some nodes to simplify the network. One possible solution is that we can search the intersections between middle nodes and pathway/function lists. Then we merge these nodes and replace them with the ID of the corresponding pathway or function.
The argument for merging is “--merge” or “-m”, 

```Python run.py -l 2  -m ```

![merge](https://github.com/bruce803/NIP/blob/master/result/merge.png)

###	filter
Suppose there are still lots of intersections after merging. And we do not plan to keep all the intersections because some of the intersections are less important. Hence, we propose two rules to filter the intersections. The first rule is the length of the intersection between pathway/function and gene list in the middle layer. We discard the intersections with length less than 3. For second, user can manually filter out some intersections by providing NIP a fdr. The NIP can calculate the p value for each intersection. Then it can filter out some less significant intersections according to the threshold of p value or FDR from user.  
For example, if user want to filter out some less important pathway or function after merging (FDR<=0.003), they can run the following command,

```Python run.py -l 3 -m  -p  0.003  ```

###	delete 
If users hope to delete one or some special nodes in the visualized network. NIP toolbox also provide the delete operation. User can just specify the ID of the nodes they want to delete from the current visualization. 
For example, suppose user want to delete function “PTK2B” and gene “SYK” from the current picture, the command is,

```Python run.py -d PTK2B SYK```

The delete operation can accept one or more arguments. So, we can delete one gene or function node, or delete some nodes at a time.

![Delete the gene “IRF4” in Figure 2 (left panel) ](https://github.com/bruce803/NIP/blob/master/result/delete.png)

