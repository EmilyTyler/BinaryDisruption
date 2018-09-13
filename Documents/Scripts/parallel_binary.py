import pycuda.driver as cuda
import pycuda.autoinit
from pycuda.compiler import SourceModule

import numpy

#Create 4x4 array of random numbers, most nvidia devices only support single precision so we have to convert (I should check if I can use double)
n_cols = 5
n_rows = 4
a = numpy.random.randn(n_rows, n_cols)
a = a.astype(numpy.float32)

#Allocate memory
a_gpu = cuda.mem_alloc(a.nbytes)

#Transfer data to the GPU
cuda.memcpy_htod(a_gpu, a)

#Write c code to perform operations on the cpu
mod = SourceModule("""
        __global__ void doublify(float *a, int lda)
        {
                int idx = threadIdx.x + threadIdx.y*lda;
                a[idx] *= 2;
        }
        """)

#Call function specifying a_gpu as the argument and a block size of 4x4
func = mod.get_function("doublify")
lda = numpy.int32(a.shape[-1])
func(a_gpu, lda, block=(n_cols,n_rows,1))

#Fetch data back from the gpu and display it alongside the original data
a_doubled = numpy.empty_like(a)
cuda.memcpy_dtoh(a_doubled, a_gpu)
print(a_doubled-2.0*a)
