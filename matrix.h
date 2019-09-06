#ifndef MATRIX_H
#define MATRIX_H

double Ralpha;

double MatrixIdentity[4][4]={{1,0,0,0},
                             {0,1,0,0},
                             {0,0,1,0},
                             {0,0,0,1}};

double MatrixTranslete[4][4]={{1,0,0,0},
                              {0,1,0,0},
                              {0,0,1,0},
                              {0,0,0,1}};

double MatrixScale[4][4]={{1,0,0,0},
                          {0,1,0,0},
                          {0,0,1,0},
                          {0,0,0,1}};

double MatrixRotatedX[4][4]={{1,0,0,0},
                             {0,1,0,0},
                             {0,0, 1,0},
                             {0,0,0,1}};

double MatrixRotatedY[4][4]={{1,0,0,0},
                             {0,1,0,0},
                             {0,0,1,0},
                             {0,0,0,1}};

double MatrixRotatedZ[4][4]={{1,0,0,0},
                             {0,1,0,0},
                             {0,0,1,0},
                             {0,0,0,1}};

void initializedMatrixTranslate(double Tx,double Ty,double Tz)
{
    MatrixTranslete[0][3]=Tx;
    MatrixTranslete[1][3]=Ty;
    MatrixTranslete[2][3]=Tz;
}

void initializedMatrixScale(double Sx,double Sy,double Sz)
{
    MatrixScale[0][0]=Sx;
    MatrixScale[1][1]=Sy;
    MatrixScale[2][2]=Sz;
}

void initializedMatrixRotatedX(double Ralpha)
{
    MatrixRotatedX[1][1]=cos(Ralpha);
    MatrixRotatedX[1][2]=-sin(Ralpha);
    MatrixRotatedX[2][1]=sin(Ralpha);
    MatrixRotatedX[2][2]=cos(Ralpha);
}

void initializedMatrixRotatedY(double Ralpha)
{
    MatrixRotatedY[0][0]=cos(Ralpha);
    MatrixRotatedY[0][2]=sin(Ralpha);
    MatrixRotatedY[2][0]=-sin(Ralpha);
    MatrixRotatedY[2][2]=cos(Ralpha);
}

void initializedMatrixRotatedZ(double Ralpha)
{
    MatrixRotatedZ[0][0]=cos(Ralpha);
    MatrixRotatedZ[0][1]=-sin(Ralpha);
    MatrixRotatedZ[1][0]=sin(Ralpha);
    MatrixRotatedZ[1][1]=cos(Ralpha);
}

void productMatrix(double (&a)[4][4], double (&b)[4][4], double (&rez)[4][4])
{
    double rezMatrix[4][4]={{0,0,0,0},
                            {0,0,0,0},
                            {0,0,0,0},
                            {0,0,0,0}};
    for(unsigned int row=0;row<4;++row)
    {
        for(unsigned int col=0;col<4;++col)
        {
            for(unsigned int inner=0;inner<4;++inner)
            {
                rezMatrix[row][col]+=a[row][inner]*b[inner][col];
            }
        }
    }
    for(unsigned int i=0;i<4;++i)
    {
        for(unsigned int j=0;j<4;++j)
        {
            rez[i][j]=rezMatrix[i][j];
        }
    }
}

#endif // MATRIX_H
