#ifndef _CFrame_H
#define _CFrame_H

    typedef void (*applyfunc)(struct _CFrame *frame, double t, double *qp, double *qp_dot);
    typedef struct _CFrame CFrame;

    struct _CFrame {

        // pointer to parameter array for the frame -- can contain things
        //      like the pattern speed, or the translation vector
        double *parameters;
        applyfunc apply;

    };
#endif
