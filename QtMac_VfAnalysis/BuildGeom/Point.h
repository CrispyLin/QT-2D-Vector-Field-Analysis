//point.h
/*
Created and Modified by Guoning Chen
*/
#pragma once
#ifndef __POINT_H__
#define __POINT_H__

#include <QDebug>

class Point2Df
{
public:
    float x, y;
    unsigned char type;

    Point2Df(float ix = 0, float iy = 0) { x = ix; y = iy;}
};


class Point2Dd{
public:
    double x, y;
    unsigned char type;

    Point2Dd(double ix = 0, double iy = 0) { x = ix; y = iy;}
};

class Point3Df{
public:
    float x, y, z;
    unsigned char type;

    Point3Df(float ix = 0, float iy = 0, float iz = 0) { x = ix; y = iy; z = iz;}
};

class Point3Dd{
public:
    double x, y, z;
    unsigned char type;

    Point3Dd(double ix = 0, double iy = 0, double iz = 0) { x = ix; y = iy; z = iz;}
};


typedef struct point3
{
    double p[3];
}point3;


class EdgeSamplePt{
public:
    int which_edge;
    double alpha_;   // the alpha value along the edge
    int end_tri;     // the end triangle containing the image of the sample
    int backward;//0=forward, 1=backward

    EdgeSamplePt(int edgeid, double alpha) { which_edge = edgeid; alpha_ = alpha; }

    void get_sample_coords(double x[3]);  // compute the actual 3D coordinates of this sample point
};




class EdgeSamplePt_List{
public:
    EdgeSamplePt **samples;
    int num;
    int curMax;

    EdgeSamplePt_List(int initsize = 1000)
    {
        if (initsize == 0)
        {
            samples = nullptr;
            num = curMax = 0;
            return;
        }

        samples = new EdgeSamplePt *[initsize];

        if (samples == nullptr)
        {
            qDebug() << "failed to reallocate memory for EdgeSamplePt!";
            exit(-1);
        }

        int i;

        for (i=0; i<initsize; i++)
            samples[i] = nullptr;
        num = 0;
        curMax = initsize;
    }

    ~EdgeSamplePt_List()
    {
        finalize();
    }

    inline void finalize()
    {
        if (samples == nullptr)
            return;

        int i;
        for (i=0; i<num; i++)
        {
            if (samples[i] == nullptr)
                continue;

            delete samples[i];
        }

        delete [] samples;
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(EdgeSamplePt *p)
    {
        if(isFull ())
            if(!extend(30000))
                return false;             //if not enough memory available, return false
        samples[num] = p;
        num++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        num --;
        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
    {
        if(num == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(num >= curMax) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 1000)
    {
        EdgeSamplePt **temp = samples;
        samples = new EdgeSamplePt * [curMax + step];
        if( samples == nullptr)
        {
            //fail
            samples = temp;
            return false;
        }

        int i;
        for(i = 0; i < curMax; i++)
            samples[i] = temp[i];
        for(i = curMax; i < curMax+step; i++)
            samples[i] = nullptr;
        curMax += step;
        delete [] temp;
        return true;
    }

    inline void reset()
    {
        num = 0;
    }


    void display();

    void display_sel_edges(int,bool,int);

    void create_display_list();

};



/*
   This is a point that store the image of a seed when doing \tau map computation
*/
typedef struct tauPt3
{
    double x = 0., y = 0., z = 0.;
    int which_tri = -1;
    double tau = 0.;
}tauPt3;

class TauParticle_List{
public:
    tauPt3 **particles = nullptr;
    int num = 0;
    int curMax;

    TauParticle_List(int initsize = 1000)
    {
        if (initsize == 0)
        {
            particles = nullptr;
            num = curMax = 0;
            return;
        }

        particles = new tauPt3 *[initsize];

        if (particles == nullptr)
        {
            qDebug() << "failed to reallocate memory for particles!";
            exit(-1);
        }

        int i;

        for (i=0; i<initsize; i++)
            particles[i] = nullptr;
        num = 0;
        curMax = initsize;
    }

    ~TauParticle_List()
    {
        finalize();
    }

    inline void finalize()
    {
        if (particles == nullptr)
            return;

        int i;
        for (i=0; i<num; i++)
        {
            if (particles[i] == nullptr)
                continue;

            delete particles[i];
        }

        delete [] particles;
    }

    //add a new vertex to the end of the list, if it succeeds, return true
    inline bool append(tauPt3 *p)
    {
        if(isFull ())
            if(!extend(30000))
                return false;             //if not enough memory available, return false
        particles[num] = p;
        num++;
        return true;
    }

    inline bool del_End() //delete the vertex at the end of the list
    {
        if(isEmpty())  return false;
        num --;
        return true;
    }

    inline bool isEmpty()  //judge whether the list is empty
    {
        if(num == 0)   return true;
        return false;
    }

    inline bool isFull()
    {
        if(num >= curMax) return true;
        return false;
    }

    //extend the original list, if it succeeds, return true
    inline bool extend(int step = 1000)
    {
        tauPt3 **temp = particles;
        particles = new tauPt3 * [curMax + step];
        if( particles == nullptr)
        {
            //fail
            particles = temp;
            return false;
        }

        int i;
        for(i = 0; i < curMax; i++)
            particles[i] = temp[i];
        for(i = curMax; i < curMax+step; i++)
            particles[i] = nullptr;
        curMax += step;
        delete [] temp;
        return true;
    }

    inline void reset()
    {
        num = 0;
    }


};

#endif /* __POINT_H__ */
