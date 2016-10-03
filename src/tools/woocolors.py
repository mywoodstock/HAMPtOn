
HEX = ['0','1','2','3','4','5','6','7','8','9','A','B','C','D','E','F']	

def hex(n):
	n = int(n)
	n = max(0, min(n, 255))
	return HEX[(n-n%16)/16] + HEX[n%16]

def rgb2hex(R, G, B):
	return '#' + hex(R) + hex(G) + hex(B)

## some simple ones  
vdgray  = '#222222'
dgray   = '#444444'
mgray   = '#888888'
lgray   = '#cccccc'
vlgray  = '#eeeeee'

## bruegel color palette
bruegel	= [ "#d6553d",
            "#755e24",
            "#514639",
            "#986b5c",
            "#a5c596",
            "#749867",
            "#bbae76",
            "#fff2c2",
            "#e8df96",
            "#909071" ]

## raphael color palette
raphael	= [ "#76a79f",
            "#e6bfa9",
            "#f3e6bd",
            "#f4e7cf",
            "#8e3649",
            "#7e8ba9",
            "#70a56c",
            "#b76669",
            "#546288",
            "#446224",
            "#bfc0c1" ]

## motherwell color palette
motherwell = [ "#928179",
               "#6a9cad",
               "#c19f7f",
               "#b1033d",
               "#8f9061",
               "#a5a677",
               "#7e6b71",
               "#7b7b6b" ]

## rich color palette
rich = [ rgb2hex(131, 16, 29),
         rgb2hex(161, 17, 53),
         rgb2hex(163, 47, 117),
         rgb2hex(126, 101, 101),
         rgb2hex(197, 107, 35),
         rgb2hex(220, 161, 60),
         rgb2hex(0, 110, 125),
         rgb2hex(0, 96, 75),
         rgb2hex(53, 58, 144),
         rgb2hex(77, 8, 103),
         rgb2hex(61, 17, 123),
         rgb2hex(41, 50, 90) ]


## cb2.0 qualitative palettes (colorbrewer 2.0)

cb2q1 = [ ## 7 colors
          rgb2hex(127, 201, 127),
          rgb2hex(190, 174, 212),
          rgb2hex(253, 192, 134),
          rgb2hex(255, 255, 153),
          rgb2hex(56, 108, 176),
          rgb2hex(240, 2, 127),
          rgb2hex(191, 91, 23),
          rgb2hex(102, 102, 102)
        ]
        
cb2q2 = [ ## 7 colors
          rgb2hex(141, 211, 199),
          rgb2hex(255, 255, 179),
          rgb2hex(190, 186, 218),
          rgb2hex(251, 128, 114),
          rgb2hex(128, 177, 211),
          rgb2hex(253, 180, 98),
          rgb2hex(179, 222, 105),
          rgb2hex(252, 205, 229)
        ]

cb2q3 = [ ## 8 colors
          rgb2hex(141,211,199),
          rgb2hex(255,255,179),
          rgb2hex(190,186,218),
          rgb2hex(251,128,114),
          rgb2hex(128,177,211),
          rgb2hex(253,180,98),
          rgb2hex(179,222,105),
          rgb2hex(252,205,229),
          rgb2hex(217,217,217),
          rgb2hex(188,128,189)
        ]