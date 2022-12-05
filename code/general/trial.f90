        program trial

        implicit none
        real const, slope, corr, x(2), y(2), x3, y3, y4
        integer j8,j9

        j8=2
        j9=0

        !-5.74389328E-31   5.17689969E-31   0.00000000       
        x(1) = 1.12786233
        x(2) = 1.12941372
        y(1) = 8.46789543E-33
        y(2) = 9.26431269E-33
        !1.12786233       8.46789543E-33   1.12941372       9.26431269E-33   

       x3 = 1.12792969
       y3 = 9.52855721E-33

       call linreg(j8,j9,const,slope,corr,x,y)

       y4 = slope*x3 + const

       print *, y3, y4



        endprogram trial
