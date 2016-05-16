# VTP_source_code

The source code of our paper:
Fast and Exact Discrete Geodesic Computation Based on Triangle-Oriented Wavefront Propagation

Yipeng Qin\*, Xiaoguang Han\*, Hongchuan Yu, Yizhou Yu, Jianjun Zhang (* Joint first authors)

ACM Transactions on Graphics (Proceedings of SIGGRAPH 2016)

USAGE:

-m [meshFile]: input model file (current code only support .obj files)

-s [sourceIndex]: the index of source

-o [outputGeodesicDistances]: file containing the geodesic distances from source to each vertex (VertIdx 0 ~ VertIdx N) 

Example: VTP.exe -m bunny.obj -s 0 -o bunny_geoDistance.txt
