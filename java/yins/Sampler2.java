/**
 * @file   Sampler2.java
 * @author Dan Spielman <spielman@math.mit.edu>
 * @date   Thurs Feb 1 2012
 *
 * @brief  for storing pairs integers, and popping off random ones
 *         
 * 
 */


package yins;

import java.lang.*;
import java.util.*;
import java.io.*;

/**
 * For storing random objects (in this case pairs of non-neg integers) 
 * and popping off random ones.  Is a matlab code by the same name.
 *
 * is now modified so that it can grow.
 *
 */
public class Sampler2
{

    public int[] array1;  // the data
    public int[] array2;  // the data

    public int last; // the last place we stored an item


    Random rand;

    public Sampler2(int initsize, long seed) {
        rand = new Random(seed);
        this.array1 = new int[initsize];
        this.array2 = new int[initsize];
        this.last = -1;
    }
    
    public Sampler2(int initsize) {
        rand = new Random(0);
        this.array1 = new int[initsize];
        this.array2 = new int[initsize];
        this.last = -1;
    }

    public void add(int item1, int item2) {

        int val1, val2, ind;
        
        this.last = this.last + 1;

        if (this.last == 0) {
            this.array1[0] = item1;
            this.array2[0] = item2;
        }
        else {
        
            ind = rand.nextInt(this.last);
        
            if (ind != this.last) {
                val1 = this.array1[ind];
                val2 = this.array2[ind];

                this.array1[ind] = item1;
                this.array2[ind] = item2;

                this.array1[this.last] = val1;
                this.array2[this.last] = val2;
            }
            else {
                this.array1[this.last] = item1;
                this.array2[this.last] = item2;
            }

            if (this.last == this.array1.length-1) {
                int[] ar1 = new int[2*this.array1.length];
                int[] ar2 = new int[2*this.array1.length];
                for (int i = 0; i <= this.last; i++) {
                    ar1[i] = this.array1[i];
                    ar2[i] = this.array2[i];
                }
                this.array1 = ar1;
                this.array2 = ar2;
            }
        }
    }

        
    public int[] poprand() {

        int[] out;

        out = new int[2];
        
        out[0] = this.array1[this.last];
        out[1] = this.array2[this.last];

        this.last = this.last - 1;

        return out;
    }

}

