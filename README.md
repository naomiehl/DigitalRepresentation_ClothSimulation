# Cloth-Simulation

Cloth simulation is one of the most difficult aspects of computer graphics. it is a deceptively simple real-world item, but in fact presents very complex internal and environmental interactions. 

The simulation of cloth requires more than the implementation of a simple dynamic system designed only for simple cloth objects. Many issues need to be resolved, concerning the accuracy and computational efficiency, robustness, stability, collision detection and response,and constraint handling 

Since the beginning, research on clothing simulation has generally relied on explicit numerical integration to advance simulation. The technique that is developed here is a revolutionary technique in clothing simulation. This article describe a cloth simulation system that is much faster than previously reported simulation systems, due to the use of implicit integration methods. 

The resulting simulation system is much faster than the previously described clothing simulation systems. The goal of this paper is to demonstrate that implicit methods for cloth overcome the performance limits inherent in explicit simulation methods. They also introduce a simple and unified treatment of damping forces, a topic that has been largely ignored until now.


# About the code

/* ========================================================================= *

    Obs. 1:
        The function addDampingForcesStretch() in the SystemSolver solver()
        method was the one that was causing errors for us just before the
        submission, so at the time we wrote the report this function was not
        well implemented yet, we managed to make it work a little later.
        So if you want to run the exact simulation that we described in the
        report, just comment the line 28 of SystemSolver.cpp, where this
        function is called.

    Obs. 2:
        Also after the submission, we added a new simulation whose difference
        to the previous one is only the position where we initialise the
        vertices, this time we initilise them a little more far from what
        could be considered as the "rest state" of the cloth. We did this
        with the objectif of generating a simulation where the cloth would
        have a bigger movement tendency. To switch between one or another
        version, just uncomment the version desired in the Cloth class
        constructor. Both versions work better using the function
        addDampingForcesStretch(), as expected.

 * ========================================================================= */

