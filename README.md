Fortran实现的有限元算法，算法单元包括：Ansys中的shell63,soild45,fluid30.
并且实现了这三者之间的流固耦合问题，由于受到求解器的限制，只可以求解有限元中的静力分析。如果想计算模态分析和震型分析，则需要将流固耦合的
刚度矩阵和质量矩阵单独送到相关特征值/奇异值分解求解包，如Lapack和Arpack中求解。
