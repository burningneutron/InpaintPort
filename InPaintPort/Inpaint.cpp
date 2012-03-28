#include "Inpaint.h"
#include <math.h>
#include <assert.h>
#include <string.h>
#include <malloc.h>
#include <stdio.h>
#include <float.h>

struct MatStub
{
	int rows, cols;
	int c;
	int step;
	unsigned char *data;
};

MatStub *createMatStub32F(int rows, int cols, int c)
{
	MatStub *mat = new MatStub;
	
	mat->rows = rows;
	mat->cols = cols;
	mat->c = c;
	mat->step = cols*sizeof(float)*c;
	mat->data = (unsigned char*)(new float[rows*cols*c]);

	return mat;
}

MatStub *createMatStub8U(int rows, int cols, int c)
{
	MatStub *mat = new MatStub;

	mat->rows = rows;
	mat->cols = cols;
	mat->c = c;
	mat->step = cols*sizeof(unsigned char)*c;
	mat->data = new unsigned char[rows*cols*c];

	return mat;
}

void setMatStub8UC1(MatStub *mat, unsigned char val)
{
	int size = mat->step * mat->rows;

	unsigned char* ptr = mat->data;
	for( int i = 0; i < size; i++ ){
		*ptr++ = val;
	}
}

void setMatStub8UC1WithMask(MatStub *mat, unsigned char val, MatStub *mask)
{
	int size = mat->cols * mat->rows;

	unsigned char *ptr = mat->data;
	unsigned char *maskPtr = mask->data;
	for( int i = 0; i < size; i++ ){
		if( *maskPtr++ != 0) *ptr = val;
		ptr++;
	}
}

void setMatStub32FC1(MatStub *mat, float val)
{
	int size = mat->cols * mat->rows;

	float* ptr = (float*)(mat->data);
	for( int i = 0; i < size; i++ ){
		*ptr++ = val;
	}
}

void setMatStub32FC1WithMask(MatStub *mat, float val, MatStub *mask)
{
	int size = mat->cols * mat->rows;

	float *ptr = (float*)(mat->data);
	unsigned char *maskPtr = mask->data;
	for( int i = 0; i < size; i++ ){
		if( *maskPtr++ != 0) *ptr = val;
		ptr++;
	}
}

void subMatStub8UC1(MatStub *src1, MatStub *src2, MatStub *dst)
{
	int size = src1->cols * src1->rows;

	unsigned char *srcPtr1 = src1->data;
	unsigned char *srcPtr2 = src2->data;
	unsigned char *dstPtr = dst->data;

	for( int i = 0; i < size; i++ ){
		*dstPtr++ = *srcPtr1++ - *srcPtr2++;
	}
}

void releaseMatStub(MatStub ** mat)
{
	delete []((*mat)->data);
	delete *mat;
	*mat = 0;
}

void copyMatStub(const MatStub *src, const MatStub *dst)
{
	assert(src->rows == src->rows);
	assert(src->cols == dst->cols);

	int size = src->step * src->rows;
	memcpy(dst->data, src->data, size);
}

struct IPPoint2D32f
{
	float x, y;
};

#undef MATSTUB_ELEM_PTR_FAST
#define MATSTUB_ELEM_PTR_FAST( mat, row, col, pix_size )  \
	((mat).data + (size_t)(mat).step*(row) + (pix_size)*(col))

#define MIN(x, y) ((x) < (y) ? (x) : (y))
#define MAX(x, y) ((x) > (y) ? (x) : (y))

typedef unsigned char uchar;

#define CAST_8U(t) (uchar)(!((t) & ~255) ? (t) : (t) > 0 ? 255 : 0)

inline float
	min4( float a, float b, float c, float d )
{
	a = MIN(a,b);
	c = MIN(c,d);
	return MIN(a,c);
}

inline int round(double value)
{
	return (int)(value + (value >= 0 ? 0.5 : -0.5));
}

#define MATSTUB_ELEM( mat, elemtype, row, col )           \
	(*(elemtype*)MATSTUB_ELEM_PTR_FAST( mat, row, col, sizeof(elemtype)))

#define MATSTUB_3COLOR_ELEM(img,type,y,x,c) MATSTUB_ELEM(img,type,y,(x)*3+(c))
#define KNOWN  0  //known outside narrow band
#define BAND   1  //narrow band (known)
#define INSIDE 2  //unknown
#define CHANGE 3  //servise

typedef struct IPHeapElem
{
	float T;
	int i,j;
	struct IPHeapElem* prev;
	struct IPHeapElem* next;
}
IPHeapElem;


class IPPriorityQueueFloat
{
protected:
	IPHeapElem *mem,*empty,*head,*tail;
	int num,in;

public:

	bool Init( const MatStub* f )
	{
		int i,j;
		for( i = num = 0; i < f->rows; i++ )
		{
			for( j = 0; j < f->cols; j++ )
				num += MATSTUB_ELEM(*f,uchar,i,j)!=0;
		}
		if (num<=0) return false;
		mem = (IPHeapElem*)malloc((num+2)*sizeof(IPHeapElem));
		if (mem==NULL) return false;

		head       = mem;
		head->i    = head->j = -1;
		head->prev = NULL;
		head->next = mem+1;
		head->T    = -FLT_MAX;
		empty      = mem+1;
		for (i=1; i<=num; i++) {
			mem[i].prev   = mem+i-1;
			mem[i].next   = mem+i+1;
			mem[i].i      = -1;
			mem[i].T      = FLT_MAX;
		}
		tail       = mem+i;
		tail->i    = tail->j = -1;
		tail->prev = mem+i-1;
		tail->next = NULL;
		tail->T    = FLT_MAX;
		return true;
	}

	bool Add(const MatStub* f) {
		int i,j;
		for (i=0; i<f->rows; i++) {
			for (j=0; j<f->cols; j++) {
				if (MATSTUB_ELEM(*f,uchar,i,j)!=0) {
					if (!Push(i,j,0)) 
						return false;
				}
			}
		}
		return true;
	}

	bool Push(int i, int j, float T) {
		IPHeapElem *tmp=empty,*add=empty;
		if (empty==tail) return false;
		while (tmp->prev->T>T) tmp = tmp->prev;
		if (tmp!=empty) {
			add->prev->next = add->next;
			add->next->prev = add->prev;
			empty = add->next;
			add->prev = tmp->prev;
			add->next = tmp;
			add->prev->next = add;
			add->next->prev = add;
		} else {
			empty = empty->next;
		}
		add->i = i;
		add->j = j;
		add->T = T;
		in++;
		//      printf("push i %3d  j %3d  T %12.4e  in %4d\n",i,j,T,in);
		return true;
	}

	bool Pop(int *i, int *j) {
		IPHeapElem *tmp=head->next;
		if (empty==tmp) return false;
		*i = tmp->i;
		*j = tmp->j;
		tmp->prev->next = tmp->next;
		tmp->next->prev = tmp->prev;
		tmp->prev = empty->prev;
		tmp->next = empty;
		tmp->prev->next = tmp;
		tmp->next->prev = tmp;
		empty = tmp;
		in--;
		//      printf("pop  i %3d  j %3d  T %12.4e  in %4d\n",tmp->i,tmp->j,tmp->T,in);
		return true;
	}

	bool Pop(int *i, int *j, float *T) {
		IPHeapElem *tmp=head->next;
		if (empty==tmp) return false;
		*i = tmp->i;
		*j = tmp->j;
		*T = tmp->T;
		tmp->prev->next = tmp->next;
		tmp->next->prev = tmp->prev;
		tmp->prev = empty->prev;
		tmp->next = empty;
		tmp->prev->next = tmp;
		tmp->next->prev = tmp;
		empty = tmp;
		in--;
		//      printf("pop  i %3d  j %3d  T %12.4e  in %4d\n",tmp->i,tmp->j,tmp->T,in);
		return true;
	}

	IPPriorityQueueFloat(void) {
		num=in=0;
		mem=empty=head=tail=NULL;
	}

	~IPPriorityQueueFloat(void)
	{
		free( mem );
	}
};

inline float VectorScalMult(IPPoint2D32f v1,IPPoint2D32f v2) {
	return v1.x*v2.x+v1.y*v2.y;
}

inline float VectorLength(IPPoint2D32f v1) {
	return v1.x*v1.x+v1.y*v1.y;
}

///////////////////////////////////////////////////////////////////////////////////////////
//HEAP::iterator Heap_Iterator;
//HEAP Heap;

float FastMarching_solve(int i1,int j1,int i2,int j2, const MatStub* f, const MatStub* t)
{
	double sol, a11, a22, m12;
	a11=MATSTUB_ELEM(*t,float,i1,j1);
	a22=MATSTUB_ELEM(*t,float,i2,j2);
	m12=MIN(a11,a22);

	if( MATSTUB_ELEM(*f,uchar,i1,j1) != INSIDE )
		if( MATSTUB_ELEM(*f,uchar,i2,j2) != INSIDE )
			if( fabs(a11-a22) >= 1.0 )
				sol = 1+m12;
			else
				sol = (a11+a22+sqrt((double)(2-(a11-a22)*(a11-a22))))*0.5;
		else
			sol = 1+a11;
	else if( MATSTUB_ELEM(*f,uchar,i2,j2) != INSIDE )
		sol = 1+a22;
	else
		sol = 1+m12;

	return (float)sol;
}

/////////////////////////////////////////////////////////////////////////////////////


static void
	icvCalcFMM(const MatStub *f, MatStub *t, IPPriorityQueueFloat *Heap, bool negate) {
		int i, j, ii = 0, jj = 0, q;
		float dist;

		while (Heap->Pop(&ii,&jj)) {

			unsigned known=(negate)?CHANGE:KNOWN;
			MATSTUB_ELEM(*f,uchar,ii,jj) = (uchar)known;

			for (q=0; q<4; q++) {
				i=0; j=0;
				if     (q==0) {i=ii-1; j=jj;}
				else if(q==1) {i=ii;   j=jj-1;}
				else if(q==2) {i=ii+1; j=jj;}
				else {i=ii;   j=jj+1;}
				if ((i<=0)||(j<=0)||(i>f->rows)||(j>f->cols)) continue;

				if (MATSTUB_ELEM(*f,uchar,i,j)==INSIDE) {
					dist = min4(FastMarching_solve(i-1,j,i,j-1,f,t),
						FastMarching_solve(i+1,j,i,j-1,f,t),
						FastMarching_solve(i-1,j,i,j+1,f,t),
						FastMarching_solve(i+1,j,i,j+1,f,t));
					MATSTUB_ELEM(*t,float,i,j) = dist;
					MATSTUB_ELEM(*f,uchar,i,j) = BAND;
					Heap->Push(i,j,dist);
				}
			}
		}

		if (negate) {
			for (i=0; i<f->rows; i++) {
				for(j=0; j<f->cols; j++) {
					if (MATSTUB_ELEM(*f,uchar,i,j) == CHANGE) {
						MATSTUB_ELEM(*f,uchar,i,j) = KNOWN;
						MATSTUB_ELEM(*t,float,i,j) = -MATSTUB_ELEM(*t,float,i,j);
					}
				}
			}
		}
}


static void
	icvTeleaInpaintFMM(const MatStub *f, MatStub *t, MatStub *out, int range, IPPriorityQueueFloat *Heap ) {
		int i = 0, j = 0, ii = 0, jj = 0, k, l, q, color = 0;
		float dist;

		if (out->c == 3) {

			while (Heap->Pop(&ii,&jj)) {

				MATSTUB_ELEM(*f,uchar,ii,jj) = KNOWN;
				for(q=0; q<4; q++) {
					if     (q==0) {i=ii-1; j=jj;}
					else if(q==1) {i=ii;   j=jj-1;}
					else if(q==2) {i=ii+1; j=jj;}
					else if(q==3) {i=ii;   j=jj+1;}
					if ((i<=1)||(j<=1)||(i>t->rows-1)||(j>t->cols-1)) continue;

					if (MATSTUB_ELEM(*f,uchar,i,j)==INSIDE) {
						dist = min4(FastMarching_solve(i-1,j,i,j-1,f,t),
							FastMarching_solve(i+1,j,i,j-1,f,t),
							FastMarching_solve(i-1,j,i,j+1,f,t),
							FastMarching_solve(i+1,j,i,j+1,f,t));
						MATSTUB_ELEM(*t,float,i,j) = dist;

						for (color=0; color<=2; color++) {
							IPPoint2D32f gradI,gradT,r;
							float Ia=0,Jx=0,Jy=0,s=1.0e-20f,w,dst,lev,dir,sat;

							if (MATSTUB_ELEM(*f,uchar,i,j+1)!=INSIDE) {
								if (MATSTUB_ELEM(*f,uchar,i,j-1)!=INSIDE) {
									gradT.x=(float)((MATSTUB_ELEM(*t,float,i,j+1)-MATSTUB_ELEM(*t,float,i,j-1)))*0.5f;
								} else {
									gradT.x=(float)((MATSTUB_ELEM(*t,float,i,j+1)-MATSTUB_ELEM(*t,float,i,j)));
								}
							} else {
								if (MATSTUB_ELEM(*f,uchar,i,j-1)!=INSIDE) {
									gradT.x=(float)((MATSTUB_ELEM(*t,float,i,j)-MATSTUB_ELEM(*t,float,i,j-1)));
								} else {
									gradT.x=0;
								}
							}
							if (MATSTUB_ELEM(*f,uchar,i+1,j)!=INSIDE) {
								if (MATSTUB_ELEM(*f,uchar,i-1,j)!=INSIDE) {
									gradT.y=(float)((MATSTUB_ELEM(*t,float,i+1,j)-MATSTUB_ELEM(*t,float,i-1,j)))*0.5f;
								} else {
									gradT.y=(float)((MATSTUB_ELEM(*t,float,i+1,j)-MATSTUB_ELEM(*t,float,i,j)));
								}
							} else {
								if (MATSTUB_ELEM(*f,uchar,i-1,j)!=INSIDE) {
									gradT.y=(float)((MATSTUB_ELEM(*t,float,i,j)-MATSTUB_ELEM(*t,float,i-1,j)));
								} else {
									gradT.y=0;
								}
							}
							for (k=i-range; k<=i+range; k++) {
								int km=k-1+(k==1),kp=k-1-(k==t->rows-2);
								for (l=j-range; l<=j+range; l++) {
									int lm=l-1+(l==1),lp=l-1-(l==t->cols-2);
									if (k>0&&l>0&&k<t->rows-1&&l<t->cols-1) {
										if ((MATSTUB_ELEM(*f,uchar,k,l)!=INSIDE)&&
											((l-j)*(l-j)+(k-i)*(k-i)<=range*range)) {
												r.y     = (float)(i-k);
												r.x     = (float)(j-l);

												dst = (float)(1./(VectorLength(r)*sqrt((double)VectorLength(r))));
												lev = (float)(1./(1+fabs(MATSTUB_ELEM(*t,float,k,l)-MATSTUB_ELEM(*t,float,i,j))));

												dir=VectorScalMult(r,gradT);
												if (fabs(dir)<=0.01) dir=0.000001f;
												w = (float)fabs(dst*lev*dir);

												if (MATSTUB_ELEM(*f,uchar,k,l+1)!=INSIDE) {
													if (MATSTUB_ELEM(*f,uchar,k,l-1)!=INSIDE) {
														gradI.x=(float)((MATSTUB_3COLOR_ELEM(*out,uchar,km,lp+1,color)-MATSTUB_3COLOR_ELEM(*out,uchar,km,lm-1,color)))*2.0f;
													} else {
														gradI.x=(float)((MATSTUB_3COLOR_ELEM(*out,uchar,km,lp+1,color)-MATSTUB_3COLOR_ELEM(*out,uchar,km,lm,color)));
													}
												} else {
													if (MATSTUB_ELEM(*f,uchar,k,l-1)!=INSIDE) {
														gradI.x=(float)((MATSTUB_3COLOR_ELEM(*out,uchar,km,lp,color)-MATSTUB_3COLOR_ELEM(*out,uchar,km,lm-1,color)));
													} else {
														gradI.x=0;
													}
												}
												if (MATSTUB_ELEM(*f,uchar,k+1,l)!=INSIDE) {
													if (MATSTUB_ELEM(*f,uchar,k-1,l)!=INSIDE) {
														gradI.y=(float)((MATSTUB_3COLOR_ELEM(*out,uchar,kp+1,lm,color)-MATSTUB_3COLOR_ELEM(*out,uchar,km-1,lm,color)))*2.0f;
													} else {
														gradI.y=(float)((MATSTUB_3COLOR_ELEM(*out,uchar,kp+1,lm,color)-MATSTUB_3COLOR_ELEM(*out,uchar,km,lm,color)));
													}
												} else {
													if (MATSTUB_ELEM(*f,uchar,k-1,l)!=INSIDE) {
														gradI.y=(float)((MATSTUB_3COLOR_ELEM(*out,uchar,kp,lm,color)-MATSTUB_3COLOR_ELEM(*out,uchar,km-1,lm,color)));
													} else {
														gradI.y=0;
													}
												}
												Ia += (float)w * (float)(MATSTUB_3COLOR_ELEM(*out,uchar,km,lm,color));
												Jx -= (float)w * (float)(gradI.x*r.x);
												Jy -= (float)w * (float)(gradI.y*r.y);
												s  += w;
										}
									}
								}
							}
							sat = (float)((Ia/s+(Jx+Jy)/(sqrt(Jx*Jx+Jy*Jy)+1.0e-20f)+0.5f));
							{
								int isat = (int)(sat);
								MATSTUB_3COLOR_ELEM(*out,uchar,i-1,j-1,color) = CAST_8U(isat);
							}
						}

						MATSTUB_ELEM(*f,uchar,i,j) = BAND;
						Heap->Push(i,j,dist);
					}
				}
			}

		} else if (out->c == 1) {

			while (Heap->Pop(&ii,&jj)) {

				MATSTUB_ELEM(*f,uchar,ii,jj) = KNOWN;
				for(q=0; q<4; q++) {
					if     (q==0) {i=ii-1; j=jj;}
					else if(q==1) {i=ii;   j=jj-1;}
					else if(q==2) {i=ii+1; j=jj;}
					else if(q==3) {i=ii;   j=jj+1;}
					if ((i<=1)||(j<=1)||(i>t->rows-1)||(j>t->cols-1)) continue;

					if (MATSTUB_ELEM(*f,uchar,i,j)==INSIDE) {
						dist = min4(FastMarching_solve(i-1,j,i,j-1,f,t),
							FastMarching_solve(i+1,j,i,j-1,f,t),
							FastMarching_solve(i-1,j,i,j+1,f,t),
							FastMarching_solve(i+1,j,i,j+1,f,t));
						MATSTUB_ELEM(*t,float,i,j) = dist;

						for (color=0; color<=0; color++) {
							IPPoint2D32f gradI,gradT,r;
							float Ia=0,Jx=0,Jy=0,s=1.0e-20f,w,dst,lev,dir,sat;

							if (MATSTUB_ELEM(*f,uchar,i,j+1)!=INSIDE) {
								if (MATSTUB_ELEM(*f,uchar,i,j-1)!=INSIDE) {
									gradT.x=(float)((MATSTUB_ELEM(*t,float,i,j+1)-MATSTUB_ELEM(*t,float,i,j-1)))*0.5f;
								} else {
									gradT.x=(float)((MATSTUB_ELEM(*t,float,i,j+1)-MATSTUB_ELEM(*t,float,i,j)));
								}
							} else {
								if (MATSTUB_ELEM(*f,uchar,i,j-1)!=INSIDE) {
									gradT.x=(float)((MATSTUB_ELEM(*t,float,i,j)-MATSTUB_ELEM(*t,float,i,j-1)));
								} else {
									gradT.x=0;
								}
							}
							if (MATSTUB_ELEM(*f,uchar,i+1,j)!=INSIDE) {
								if (MATSTUB_ELEM(*f,uchar,i-1,j)!=INSIDE) {
									gradT.y=(float)((MATSTUB_ELEM(*t,float,i+1,j)-MATSTUB_ELEM(*t,float,i-1,j)))*0.5f;
								} else {
									gradT.y=(float)((MATSTUB_ELEM(*t,float,i+1,j)-MATSTUB_ELEM(*t,float,i,j)));
								}
							} else {
								if (MATSTUB_ELEM(*f,uchar,i-1,j)!=INSIDE) {
									gradT.y=(float)((MATSTUB_ELEM(*t,float,i,j)-MATSTUB_ELEM(*t,float,i-1,j)));
								} else {
									gradT.y=0;
								}
							}
							for (k=i-range; k<=i+range; k++) {
								int km=k-1+(k==1),kp=k-1-(k==t->rows-2);
								for (l=j-range; l<=j+range; l++) {
									int lm=l-1+(l==1),lp=l-1-(l==t->cols-2);
									if (k>0&&l>0&&k<t->rows-1&&l<t->cols-1) {
										if ((MATSTUB_ELEM(*f,uchar,k,l)!=INSIDE)&&
											((l-j)*(l-j)+(k-i)*(k-i)<=range*range)) {
												r.y     = (float)(i-k);
												r.x     = (float)(j-l);

												dst = (float)(1./(VectorLength(r)*sqrt(VectorLength(r))));
												lev = (float)(1./(1+fabs(MATSTUB_ELEM(*t,float,k,l)-MATSTUB_ELEM(*t,float,i,j))));

												dir=VectorScalMult(r,gradT);
												if (fabs(dir)<=0.01) dir=0.000001f;
												w = (float)fabs(dst*lev*dir);

												if (MATSTUB_ELEM(*f,uchar,k,l+1)!=INSIDE) {
													if (MATSTUB_ELEM(*f,uchar,k,l-1)!=INSIDE) {
														gradI.x=(float)((MATSTUB_ELEM(*out,uchar,km,lp+1)-MATSTUB_ELEM(*out,uchar,km,lm-1)))*2.0f;
													} else {
														gradI.x=(float)((MATSTUB_ELEM(*out,uchar,km,lp+1)-MATSTUB_ELEM(*out,uchar,km,lm)));
													}
												} else {
													if (MATSTUB_ELEM(*f,uchar,k,l-1)!=INSIDE) {
														gradI.x=(float)((MATSTUB_ELEM(*out,uchar,km,lp)-MATSTUB_ELEM(*out,uchar,km,lm-1)));
													} else {
														gradI.x=0;
													}
												}
												if (MATSTUB_ELEM(*f,uchar,k+1,l)!=INSIDE) {
													if (MATSTUB_ELEM(*f,uchar,k-1,l)!=INSIDE) {
														gradI.y=(float)((MATSTUB_ELEM(*out,uchar,kp+1,lm)-MATSTUB_ELEM(*out,uchar,km-1,lm)))*2.0f;
													} else {
														gradI.y=(float)((MATSTUB_ELEM(*out,uchar,kp+1,lm)-MATSTUB_ELEM(*out,uchar,km,lm)));
													}
												} else {
													if (MATSTUB_ELEM(*f,uchar,k-1,l)!=INSIDE) {
														gradI.y=(float)((MATSTUB_ELEM(*out,uchar,kp,lm)-MATSTUB_ELEM(*out,uchar,km-1,lm)));
													} else {
														gradI.y=0;
													}
												}
												Ia += (float)w * (float)(MATSTUB_ELEM(*out,uchar,km,lm));
												Jx -= (float)w * (float)(gradI.x*r.x);
												Jy -= (float)w * (float)(gradI.y*r.y);
												s  += w;
										}
									}
								}
							}
							sat = (float)((Ia/s+(Jx+Jy)/(sqrt(Jx*Jx+Jy*Jy)+1.0e-20f)+0.5f));
							{
								int isat = round(sat);
								MATSTUB_ELEM(*out,uchar,i-1,j-1) = CAST_8U(isat);
							}
						}

						MATSTUB_ELEM(*f,uchar,i,j) = BAND;
						Heap->Push(i,j,dist);
					}
				}
			}
		}
}


static void
	icvNSInpaintFMM(const MatStub *f, MatStub *t, MatStub *out, int range, IPPriorityQueueFloat *Heap) {
		int i = 0, j = 0, ii = 0, jj = 0, k, l, q, color = 0;
		float dist;

		if (out->c == 3) {

			while (Heap->Pop(&ii,&jj)) {

				MATSTUB_ELEM(*f,uchar,ii,jj) = KNOWN;
				for(q=0; q<4; q++) {
					if     (q==0) {i=ii-1; j=jj;}
					else if(q==1) {i=ii;   j=jj-1;}
					else if(q==2) {i=ii+1; j=jj;}
					else if(q==3) {i=ii;   j=jj+1;}
					if ((i<=1)||(j<=1)||(i>t->rows-1)||(j>t->cols-1)) continue;

					if (MATSTUB_ELEM(*f,uchar,i,j)==INSIDE) {
						dist = min4(FastMarching_solve(i-1,j,i,j-1,f,t),
							FastMarching_solve(i+1,j,i,j-1,f,t),
							FastMarching_solve(i-1,j,i,j+1,f,t),
							FastMarching_solve(i+1,j,i,j+1,f,t));
						MATSTUB_ELEM(*t,float,i,j) = dist;

						for (color=0; color<=2; color++) {
							IPPoint2D32f gradI,r;
							float Ia=0,s=1.0e-20f,w,dst,dir;

							for (k=i-range; k<=i+range; k++) {
								int km=k-1+(k==1),kp=k-1-(k==f->rows-2);
								for (l=j-range; l<=j+range; l++) {
									int lm=l-1+(l==1),lp=l-1-(l==f->cols-2);
									if (k>0&&l>0&&k<f->rows-1&&l<f->cols-1) {
										if ((MATSTUB_ELEM(*f,uchar,k,l)!=INSIDE)&&
											((l-j)*(l-j)+(k-i)*(k-i)<=range*range)) {
												r.y=(float)(k-i);
												r.x=(float)(l-j);

												dst = 1/(VectorLength(r)*VectorLength(r)+1);

												if (MATSTUB_ELEM(*f,uchar,k+1,l)!=INSIDE) {
													if (MATSTUB_ELEM(*f,uchar,k-1,l)!=INSIDE) {
														gradI.x=(float)(abs(MATSTUB_3COLOR_ELEM(*out,uchar,kp+1,lm,color)-MATSTUB_3COLOR_ELEM(*out,uchar,kp,lm,color))+
															abs(MATSTUB_3COLOR_ELEM(*out,uchar,kp,lm,color)-MATSTUB_3COLOR_ELEM(*out,uchar,km-1,lm,color)));
													} else {
														gradI.x=(float)(abs(MATSTUB_3COLOR_ELEM(*out,uchar,kp+1,lm,color)-MATSTUB_3COLOR_ELEM(*out,uchar,kp,lm,color)))*2.0f;
													}
												} else {
													if (MATSTUB_ELEM(*f,uchar,k-1,l)!=INSIDE) {
														gradI.x=(float)(abs(MATSTUB_3COLOR_ELEM(*out,uchar,kp,lm,color)-MATSTUB_3COLOR_ELEM(*out,uchar,km-1,lm,color)))*2.0f;
													} else {
														gradI.x=0;
													}
												}
												if (MATSTUB_ELEM(*f,uchar,k,l+1)!=INSIDE) {
													if (MATSTUB_ELEM(*f,uchar,k,l-1)!=INSIDE) {
														gradI.y=(float)(abs(MATSTUB_3COLOR_ELEM(*out,uchar,km,lp+1,color)-MATSTUB_3COLOR_ELEM(*out,uchar,km,lm,color))+
															abs(MATSTUB_3COLOR_ELEM(*out,uchar,km,lm,color)-MATSTUB_3COLOR_ELEM(*out,uchar,km,lm-1,color)));
													} else {
														gradI.y=(float)(abs(MATSTUB_3COLOR_ELEM(*out,uchar,km,lp+1,color)-MATSTUB_3COLOR_ELEM(*out,uchar,km,lm,color)))*2.0f;
													}
												} else {
													if (MATSTUB_ELEM(*f,uchar,k,l-1)!=INSIDE) {
														gradI.y=(float)(abs(MATSTUB_3COLOR_ELEM(*out,uchar,km,lm,color)-MATSTUB_3COLOR_ELEM(*out,uchar,km,lm-1,color)))*2.0f;
													} else {
														gradI.y=0;
													}
												}

												gradI.x=-gradI.x;
												dir=VectorScalMult(r,gradI);

												if (fabs(dir)<=0.01) {
													dir=0.000001f;
												} else {
													dir = (float)fabs(VectorScalMult(r,gradI)/sqrt(VectorLength(r)*VectorLength(gradI)));
												}
												w = dst*dir;
												Ia += (float)w * (float)(MATSTUB_3COLOR_ELEM(*out,uchar,km,lm,color));
												s  += w;
										}
									}
								}
							}
							{
								int out_val = round((double)Ia/s);
								MATSTUB_3COLOR_ELEM(*out,uchar,i-1,j-1,color) = CAST_8U(out_val);
							}
						}

						MATSTUB_ELEM(*f,uchar,i,j) = BAND;
						Heap->Push(i,j,dist);
					}
				}
			}

		} else if (out->c == 3) {

			while (Heap->Pop(&ii,&jj)) {

				MATSTUB_ELEM(*f,uchar,ii,jj) = KNOWN;
				for(q=0; q<4; q++) {
					if     (q==0) {i=ii-1; j=jj;}
					else if(q==1) {i=ii;   j=jj-1;}
					else if(q==2) {i=ii+1; j=jj;}
					else if(q==3) {i=ii;   j=jj+1;}
					if ((i<=1)||(j<=1)||(i>t->rows-1)||(j>t->cols-1)) continue;

					if (MATSTUB_ELEM(*f,uchar,i,j)==INSIDE) {
						dist = min4(FastMarching_solve(i-1,j,i,j-1,f,t),
							FastMarching_solve(i+1,j,i,j-1,f,t),
							FastMarching_solve(i-1,j,i,j+1,f,t),
							FastMarching_solve(i+1,j,i,j+1,f,t));
						MATSTUB_ELEM(*t,float,i,j) = dist;

						{
							IPPoint2D32f gradI,r;
							float Ia=0,s=1.0e-20f,w,dst,dir;

							for (k=i-range; k<=i+range; k++) {
								int km=k-1+(k==1),kp=k-1-(k==t->rows-2);
								for (l=j-range; l<=j+range; l++) {
									int lm=l-1+(l==1),lp=l-1-(l==t->cols-2);
									if (k>0&&l>0&&k<t->rows-1&&l<t->cols-1) {
										if ((MATSTUB_ELEM(*f,uchar,k,l)!=INSIDE)&&
											((l-j)*(l-j)+(k-i)*(k-i)<=range*range)) {
												r.y=(float)(i-k);
												r.x=(float)(j-l);

												dst = 1/(VectorLength(r)*VectorLength(r)+1);

												if (MATSTUB_ELEM(*f,uchar,k+1,l)!=INSIDE) {
													if (MATSTUB_ELEM(*f,uchar,k-1,l)!=INSIDE) {
														gradI.x=(float)(abs(MATSTUB_ELEM(*out,uchar,kp+1,lm)-MATSTUB_ELEM(*out,uchar,kp,lm))+
															abs(MATSTUB_ELEM(*out,uchar,kp,lm)-MATSTUB_ELEM(*out,uchar,km-1,lm)));
													} else {
														gradI.x=(float)(abs(MATSTUB_ELEM(*out,uchar,kp+1,lm)-MATSTUB_ELEM(*out,uchar,kp,lm)))*2.0f;
													}
												} else {
													if (MATSTUB_ELEM(*f,uchar,k-1,l)!=INSIDE) {
														gradI.x=(float)(abs(MATSTUB_ELEM(*out,uchar,kp,lm)-MATSTUB_ELEM(*out,uchar,km-1,lm)))*2.0f;
													} else {
														gradI.x=0;
													}
												}
												if (MATSTUB_ELEM(*f,uchar,k,l+1)!=INSIDE) {
													if (MATSTUB_ELEM(*f,uchar,k,l-1)!=INSIDE) {
														gradI.y=(float)(abs(MATSTUB_ELEM(*out,uchar,km,lp+1)-MATSTUB_ELEM(*out,uchar,km,lm))+
															abs(MATSTUB_ELEM(*out,uchar,km,lm)-MATSTUB_ELEM(*out,uchar,km,lm-1)));
													} else {
														gradI.y=(float)(abs(MATSTUB_ELEM(*out,uchar,km,lp+1)-MATSTUB_ELEM(*out,uchar,km,lm)))*2.0f;
													}
												} else {
													if (MATSTUB_ELEM(*f,uchar,k,l-1)!=INSIDE) {
														gradI.y=(float)(abs(MATSTUB_ELEM(*out,uchar,km,lm)-MATSTUB_ELEM(*out,uchar,km,lm-1)))*2.0f;
													} else {
														gradI.y=0;
													}
												}

												gradI.x=-gradI.x;
												dir=VectorScalMult(r,gradI);

												if (fabs(dir)<=0.01) {
													dir=0.000001f;
												} else {
													dir = (float)fabs(VectorScalMult(r,gradI)/sqrt(VectorLength(r)*VectorLength(gradI)));
												}
												w = dst*dir;
												Ia += (float)w * (float)(MATSTUB_ELEM(*out,uchar,km,lm));
												s  += w;
										}
									}
								}
							}
							{
								int out_val = round((double)Ia/s);
								MATSTUB_ELEM(*out,uchar,i-1,j-1) = CAST_8U(out_val);
							}
						}

						MATSTUB_ELEM(*f,uchar,i,j) = BAND;
						Heap->Push(i,j,dist);
					}
				}
			}

		}
}

#define SET_BORDER1_C1(image,type,value) {\
	int i,j;\
	for(j=0; j<image->cols; j++) {\
	MATSTUB_ELEM(*image,type,0,j) = value;\
	}\
	for (i=1; i<image->rows-1; i++) {\
	MATSTUB_ELEM(*image,type,i,0) = MATSTUB_ELEM(*image,type,i,image->cols-1) = value;\
	}\
	for(j=0; j<image->cols; j++) {\
	MATSTUB_ELEM(*image,type,erows-1,j) = value;\
	}\
}

#define COPY_MASK_BORDER1_C1(src,dst,type) {\
	int i,j;\
	for (i=0; i<src->rows; i++) {\
	for(j=0; j<src->cols; j++) {\
	if (MATSTUB_ELEM(*src,type,i,j)!=0)\
	MATSTUB_ELEM(*dst,type,i+1,j+1) = INSIDE;\
	}\
	}\
}


void dilate(MatStub *src, MatStub *dst)
{
	int width = src->cols;
	int height = src->rows;

	unsigned char *srcPtr = src->data + width + 1;
	unsigned char *dstPtr = dst->data + width + 1;
	unsigned char *p1, *p2, *p3, *p4, *p6, *p7, *p8, *p9;

	p1=src->data;p2=p1+1;p3=p2+1;
	p4=p1+width; p6=p4+2;
	p7=p4+width; p8=p7+1; p9=p8+1;

	for (int i = 1; i < height - 1; ++i) {
		for (int j = 1; j < width - 1; ++j) {
			if (*p1 != 0
				|| *p2 != 0
				|| *p3 != 0
				|| *p4 != 0
				|| *p6 != 0
				|| *p7 != 0
				|| *p8 != 0
				|| *p9 != 0)
			{
				 *dstPtr = INSIDE;
			}
			++p1;++p2;++p3;++p4;++srcPtr;++p6;++p7;++p8;++p9;
			++dstPtr;
		}
		p1+=2;p2+=2;p3+=2;p4+=2;srcPtr+=2;p6+=2;p7+=2;p8+=2;p9+=2;
		dstPtr+=2;
	}
}

void inpaint( MatStub *input_img, MatStub *inpaint_mask, MatStub *output_img, double inpaintRange)
{
	MatStub *t, *mask, *band;
	IPPriorityQueueFloat* Heap;
	int range=round(inpaintRange);    
	int erows, ecols;

	range = MAX(range,1);
	range = MIN(range,100);

	ecols = input_img->cols + 2;
	erows = input_img->rows + 2;

	t = createMatStub32F(erows, ecols, 1);
	band = createMatStub8U(erows, ecols, 1);
	mask = createMatStub8U(erows, ecols, 1);

	copyMatStub( input_img, output_img );
	setMatStub8UC1(mask,KNOWN);
	COPY_MASK_BORDER1_C1(inpaint_mask,mask,uchar);
	SET_BORDER1_C1(mask,uchar,0);
 	setMatStub32FC1(t, 1.0e6f);

	//cvDilate(mask,band,el_cross,1);   // image with narrow band
	setMatStub8UC1(band, 0);
	dilate(mask, band);
	//
	Heap=new IPPriorityQueueFloat;
	if (!Heap->Init(band)){
		assert(0);
		return;
	}
	subMatStub8UC1(band,mask,band);
	SET_BORDER1_C1(band,uchar,0);

	if (!Heap->Add(band)){
		assert(0);
		return;
	}	
	
	setMatStub32FC1WithMask(t,0,band);
	
	icvNSInpaintFMM(mask,t,output_img,range,Heap);

	releaseMatStub(&t);
	releaseMatStub(&band);
	releaseMatStub(&mask);

	delete Heap;
}

void Inpaint( unsigned char *src, unsigned char *mask, int rows, int cols, double inpaintRange)
{
	MatStub *outMat = createMatStub8U(rows, cols, 3);
	MatStub srcMat;
	srcMat.rows = rows;
	srcMat.cols = cols;
	srcMat.c = 3;
	srcMat.step = 3*cols;
	srcMat.data = src;

	MatStub maskMat;
	maskMat.rows = rows;
	maskMat.cols = cols;
	maskMat.c = 1;
	maskMat.step = cols;
	maskMat.data = mask;

	
	inpaint(&srcMat, &maskMat, outMat, inpaintRange);

	copyMatStub(outMat, &srcMat);

	releaseMatStub(&outMat);
}

void abgr2bgr(unsigned char *abgr, unsigned char *bgr, int w, int h)
{
	for( int i = 0; i < w*h; i++ ){
		abgr++;
		*bgr++ = *abgr++;
		*bgr++ = *abgr++;
		*bgr++ = *abgr++;
	}
}

void bgr2abgr(unsigned char *bgr, unsigned char *abgr, int w, int h)
{
	for( int i = 0; i < w*h; i++ ){
		abgr++;
		*abgr++ = *bgr++;
		*abgr++ = *bgr++;
		*abgr++ = *bgr++;
	}
}

void InpaintWithAlpha( unsigned char *src, unsigned char *mask, int rows, int cols, double inpaintRange)
{
	MatStub *outMat = createMatStub8U(rows, cols, 3);
	MatStub *srcMat = createMatStub8U(rows, cols, 3);
	abgr2bgr(src, srcMat->data, cols, rows);
	

	MatStub maskMat;
	maskMat.rows = rows;
	maskMat.cols = cols;
	maskMat.c = 1;
	maskMat.step = cols;
	maskMat.data = mask;


	inpaint(srcMat, &maskMat, outMat, inpaintRange);

	bgr2abgr(outMat->data, src, cols, rows);

	releaseMatStub(&srcMat);
	releaseMatStub(&outMat);
}

//void
//	inpaint( const CvArr* _input_img, const CvArr* _inpaint_mask, CvArr* _output_img,
//	double inpaintRange, int flags )
//{
//	cv::Ptr<MatStub> mask, band, f, t, out;
//	cv::Ptr<IPPriorityQueueFloat> Heap, Out;
//	cv::Ptr<IplConvKernel> el_cross, el_range;
//
//	MatStub input_hdr, mask_hdr, output_hdr;
//	MatStub* input_img, *inpaint_mask, *output_img;
//	int range=cvRound(inpaintRange);    
//	int erows, ecols;
//
//	input_img = cvGetMat( _input_img, &input_hdr );
//	inpaint_mask = cvGetMat( _inpaint_mask, &mask_hdr );
//	output_img = cvGetMat( _output_img, &output_hdr );
//
//	if( !CV_ARE_SIZES_EQ(input_img,output_img) || !CV_ARE_SIZES_EQ(input_img,inpaint_mask))
//		CV_Error( CV_StsUnmatchedSizes, "All the input and output images must have the same size" );
//
//	if( (CV_MAT_TYPE(input_img->type) != CV_8UC1 &&
//		CV_MAT_TYPE(input_img->type) != CV_8UC3) ||
//		!CV_ARE_TYPES_EQ(input_img,output_img) )
//		CV_Error( CV_StsUnsupportedFormat,
//		"Only 8-bit 1-channel and 3-channel input/output images are supported" );
//
//	if( CV_MAT_TYPE(inpaint_mask->type) != CV_8UC1 )
//		CV_Error( CV_StsUnsupportedFormat, "The mask must be 8-bit 1-channel image" );
//
//	range = MAX(range,1);
//	range = MIN(range,100);
//
//	ecols = input_img->cols + 2;
//	erows = input_img->rows + 2;
//
//	f = cvCreateMat(erows, ecols, CV_8UC1);
//	t = cvCreateMat(erows, ecols, CV_32FC1);
//	band = cvCreateMat(erows, ecols, CV_8UC1);
//	mask = cvCreateMat(erows, ecols, CV_8UC1);
//	el_cross = cvCreateStructuringElementEx(3,3,1,1,CV_SHAPE_CROSS,NULL);
//
//	cvCopy( input_img, output_img );
//	cvSet(mask,cvScalar(KNOWN,0,0,0));
//	COPY_MASK_BORDER1_C1(inpaint_mask,mask,uchar);
//	SET_BORDER1_C1(mask,uchar,0);
//	cvSet(f,cvScalar(KNOWN,0,0,0));
//	cvSet(t,cvScalar(1.0e6f,0,0,0));
//	cvDilate(mask,band,el_cross,1);   // image with narrow band
//	Heap=new IPPriorityQueueFloat;
//	if (!Heap->Init(band))
//		return;
//	cvSub(band,mask,band,NULL);
//	SET_BORDER1_C1(band,uchar,0);
//	if (!Heap->Add(band))
//		return;
//	cvSet(f,cvScalar(BAND,0,0,0),band);
//	cvSet(f,cvScalar(INSIDE,0,0,0),mask);
//	cvSet(t,cvScalar(0,0,0,0),band);
//
//	if( flags == CV_INPAINT_TELEA )
//	{
//		out = cvCreateMat(erows, ecols, CV_8UC1);
//		el_range = cvCreateStructuringElementEx(2*range+1,2*range+1,
//			range,range,CV_SHAPE_RECT,NULL);
//		cvDilate(mask,out,el_range,1);
//		cvSub(out,mask,out,NULL);
//		Out=new IPPriorityQueueFloat;
//		if (!Out->Init(out))
//			return;
//		if (!Out->Add(band))
//			return;
//		cvSub(out,band,out,NULL);
//		SET_BORDER1_C1(out,uchar,0);
//		icvCalcFMM(out,t,Out,true);
//		icvTeleaInpaintFMM(mask,t,output_img,range,Heap);
//	}
//	else if (flags == CV_INPAINT_NS) {
//		icvNSInpaintFMM(mask,t,output_img,range,Heap);
//	} else {
//		CV_Error( CV_StsBadArg, "The flags argument must be one of CV_INPAINT_TELEA or CV_INPAINT_NS" );
//	}
//}