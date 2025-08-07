# MagCoords.jl

This is a Julia library for doing coordinate transformations for Earth-centered magnetic coordinates, mainly useful for magnetosphere and ionosphere physicists, I would suppose.

It's purposefully simple, but the downside is that it has no special safety checks of the kind that might help prevent you from making mistakes; it's basically just the core math and astrometry formulas needed to do coordinate transformations. It does have tests, so you can verify that it works correctly, and includes some examples to make sure you're using it correctly, but you should do some sanity checks on your own results, too.

I made use of unicode variable names, because I felt like it made the code clearer to read. However, those are only used internally, not for the user code. The Julia syntax made what "should" be very complex code in most languages look surprisingly compact and simple (in my opinion).

See MagCoordsExamples.jl to learn how to use it.

If you have better versions of any of the astronomical functions used, please don't hesitate to suggest a change. I did my best to try to find the most accurate ones that would work as far into the future as possible, but astrometry buffs are welcome to point out better ones.

This was my first code developed in Julia, a weekend learning project (more than one weekend). Since it was not developed as part of my job, or funded by anyone but myself, I'm releasing it under the MIT License.
