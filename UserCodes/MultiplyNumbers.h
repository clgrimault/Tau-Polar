#ifndef _MULTIPLYNUMBERS_H
#define _MULTIPYNUMBERS_H

class MultiplyNumbers
{
        private:
        int _a;
        int _b;

        public:
        MultiplyNumbers ();
        ~MultiplyNumbers ();

        void setA (int a);
        void setB (int b);

        int getA () const;
        int getB () const;

        int getProduct () const;

}; // MultiplyNumbers

#endif // _ADDNUMBERS_H


