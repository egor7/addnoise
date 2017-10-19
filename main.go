package main

import (
	"fmt"
	"image"
	"log"
	"math"
	"math/cmplx"
	"math/rand"
	"os"

	"github.com/gonum/matrix/mat64"
	"github.com/mjibson/go-dsp/fft"

	"image/color"
	"image/png"
)

// "image/color"
// _ "image/gif"
// _ "image/jpeg"
// _ "image/png"

func main() {

	// s
	f_s, err := os.Open("s.png")
	if err != nil {
		log.Fatal(err)
	}
	defer f_s.Close()
	im_s, _, err := image.Decode(f_s)
	if err != nil {
		log.Fatal(err)
	}
	bounds := im_s.Bounds()
	width, height := bounds.Max.X-bounds.Min.X, bounds.Max.Y-bounds.Min.Y
	h_s := make([]float64, width*height)
	hi_s := make([][]complex128, width)
	for c := 0; c < width; c++ {
		hi_s[c] = make([]complex128, height)
		for r := 0; r < height; r++ {
			h_s[c*height+r] = float64(color.GrayModel.Convert(im_s.At(c, r)).(color.Gray).Y)
			hi_s[c][r] = complex(h_s[c*height+r], 0)
		}
	}
	fmt.Printf("read %s, size [%dx%d]\n", f_s.Name(), width, height)

	fft_s := fft.FFT2(hi_s)

	// add some noise
	rnd := rand.New(rand.NewSource(99))

	for num := 0; num < 900; num++ {
		c := rnd.Intn(width)
		r := rnd.Intn(height)

		cc := (c+width/2)%width - width/2
		rr := (r+height/2)%height - height/2

		if cc*cc+rr*rr > 10000 {
			fft_s[c][r] += complex(float64(rnd.Int63n(400000)), 0)
		}

	}

	h_c := fft.IFFT2(fft_s)

	// save c
	im_c := image.NewGray(image.Rect(0, 0, width, height))
	for c := 0; c < width; c++ {
		for r := 0; r < height; r++ {
			height := cmplx.Abs(h_c[c][r])
			if height > 255 {
				height = 255
			}
			im_c.Set(c, r, color.Gray{uint8(height)})
		}
	}
	f_c, err := os.Create("c.png")
	defer f_c.Close()
	png.Encode(f_c, im_c)
	fmt.Printf("saved %s\n", f_c.Name())
}

// side = 1 => exterior
// side = -1 => interior
//
// img in [0..width/2, 0..height/2]
func approx(img [][]complex128, width int, height int, dr int, side int) (xa, xb, xc, delta float64) {
	xa = 0
	xb = 0
	xc = 0
	delta = 0
	// (ax^2+by^2+cxy+d-f(x,y))*x^2 = 0
	//                          y^2 = 0
	//                          xy  = 0
	//                          1   = 0
	//
	// x4   x2y2 x3y  x2  | fx2
	// x2y2 y4   xy3  y2  | fy2
	// x3y  xy3  x2y2 xy  | fxy
	// x2   y2   xy   1   | f

	// instead use a plate:

	// (ax + by + c - f(x,y))^2 ~> min
	//
	// (ax + by + c - f(x,y))x
	// (ax + by + c - f(x,y))y
	// (ax + by + c - f(x,y))
	//
	// x2 xy x1 | fx
	// xy y2 y1 | fy
	// x1 y1  1 | f

	n := 0
	for c := 0; c < width/2; c++ {
		for r := 0; r < height/2; r++ {
			if side*(c*c+r*r) > side*(dr*dr) {
				n = n + 1
			}
		}
	}

	var (
		x2, y2, xy, x1, y1, I float64
		F, Fx, Fy             float64
	)
	dd := 1.0 / float64(n)
	for c := 0; c < width/2; c++ {
		for r := 0; r < height/2; r++ {
			if side*(c*c+r*r) > side*(dr*dr) {
				f_im := cmplx.Abs(img[c][r])

				x := c
				y := r

				x2 += float64(x*x) * dd
				y2 += float64(y*y) * dd
				xy += float64(x*y) * dd
				x1 += float64(x) * dd
				y1 += float64(y) * dd
				I += dd

				F += f_im * dd
				Fx += float64(x) * f_im * dd
				Fy += float64(y) * f_im * dd
			}
		}
	}

	a := mat64.NewDense(3, 3, []float64{
		x2, xy, x1,
		xy, y2, y1,
		x1, y1, I,
	})
	b := mat64.NewVector(3, []float64{Fx, Fy, F})

	var x mat64.Vector
	if err := x.SolveVec(a, b); err != nil {
		//fmt.Println("Matrix is near singular: ", err)
	}
	//fmt.Printf("x = %0.4v\n", mat64.Formatted(&x, mat64.Prefix("    ")))

	xa = x.At(0, 0)
	xb = x.At(1, 0)
	xc = x.At(2, 0)
	for c := 0; c < width/2; c++ {
		for r := 0; r < height/2; r++ {
			if side*(c*c+r*r) > side*(dr*dr) {
				f_im := cmplx.Abs(img[c][r])

				x := float64(c)
				y := float64(r)

				delta += math.Sqrt((xa*x+xb*y+xc-f_im)*(xa*x+xb*y+xc-f_im)) * dd
			}
		}
	}

	return
}

// img in [0..width/2, 0..height/2]
func abs(img [][]complex128, width, height int) ([]float64, float64) {
	h := make([]float64, width*height)
	h_max := 0.0
	for c := 0; c < width; c++ {
		for r := 0; r < height; r++ {
			cc := (c + width/2) % width
			rr := (r + height/2) % height

			h[cc*height+rr] = cmplx.Abs(img[c][r])
			if h[cc*height+rr] > h_max {
				h_max = h[cc*height+rr]
			}
		}
	}
	return h, h_max
}

// func dump(img []float64, width, height int, nm string) {
// 	f, err := os.Create(nm + ".pgm")
// 	if err != nil {
// 		log.Fatal(err)
// 	}
// 	defer f.Close()
// 	png.Encode(f_un_c, im_un_c)
// 	fmt.Printf("saved %s\n", f_un_c.Name())
// }
