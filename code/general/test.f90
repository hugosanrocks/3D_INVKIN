PROGRAM HELLO

INTEGER threads, id

threads = omp_get_num_threads()
!$OMP PARALLEL
id = omp_get_thread_num()
PRINT *, 'NUM THREADS:', THREADS
PRINT *, 'hello from thread:',id,' out of', threads
!$OMP END PARALLEL

STOP
END
