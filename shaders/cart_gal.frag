uniform vec2 iResolution;
uniform vec2 iMouse;
uniform float iTime;

//// ------ starts

// Star Nest by Pablo Roman Andrioli

// This content is under the MIT License.

#define iterations 17
#define formuparam 0.53

#define volsteps 20
#define stepsize 0.1

#define zoom 0.500
#define tile 0.850
#define speed 0.001

#define brightness 0.0015
#define darkmatter 0.300
#define distfading 0.730
#define saturation 0.850

vec4 starField(vec2 fragCoord) {
    //get coords and direction
    vec2 uv = fragCoord.xy / iResolution.xy - .5;
    uv.y *= iResolution.y / iResolution.x;
    vec3 dir = vec3(uv * zoom, 1.);
    float time = iTime * speed + .25;

    //mouse rotation
    float a1 = .5 + iMouse.x / iResolution.x * 2.;
    float a2 = .8 + iMouse.y / iResolution.y * 2.;
    mat2 rot1 = mat2(cos(a1), sin(a1), -sin(a1), cos(a1));
    mat2 rot2 = mat2(cos(a2), sin(a2), -sin(a2), cos(a2));
    dir.xz *= rot1;
    dir.xy *= rot2;
    vec3 from = vec3(1., .5, 0.5);
    from += vec3(time * 2., time, -2.);
    from.xz *= rot1;
    from.xy *= rot2;

    //volumetric rendering
    float s = 0.1, fade = 1.;
    vec3 v = vec3(0.);
    for (int r = 0; r < volsteps; r++) {
        vec3 p = from + s * dir * .5;
        p = abs(vec3(tile) - mod(p, vec3(tile * 2.))); // tiling fold
        float pa, a = pa = 0.;
        for (int i = 0; i < iterations; i++) {
            p = abs(p) / dot(p, p) - formuparam; // the magic formula
            a += abs(length(p) - pa); // absolute sum of average change
            pa = length(p);
        }
        float dm = max(0., darkmatter - a * a * .001); //dark matter
        a *= a * a; // add contrast
        if (r > 6) fade *= 1. - dm; // dark matter, don't render near
        //v+=vec3(dm,dm*.5,0.);
        v += fade;
        v += vec3(s, s * s, s * s * s * s) * a * brightness * fade; // coloring based on distance
        fade *= distfading; // distance fading
        s += stepsize;
    }
    v = mix(vec3(length(v)), v, saturation); //color adjust
    return vec4(v * .01, 1.);

}

    //// ------ galaxy

    // Created by sebastien durand - 08/2016
    //-------------------------------------------------------------------------------------
    // Based on "Dusty nebula 4" by Duke (https://www.shadertoy.com/view/MsVXWW)
    // Sliders from IcePrimitives by Bers (https://www.shadertoy.com/view/MscXzn)
    // License: Creative Commons Attribution-NonCommercial-ShareAlike 3.0 Unported License
    //-------------------------------------------------------------------------------------

    #define R(p, a) p = cos(a) * p + sin(a) * vec2(p.y, -p.x)
    #define pi 3.14159265

const vec4
colCenter = vec4(1.2, 1.5, 1.5, .25),
colEdge = vec4(.1, .1, .2, .5),
colEdge2 = vec4(.7, .54, .3, .23),
colEdge3 = vec4(.6, 1., 1.3, .25);

const float time = 10.;
vec4 sliderVal = vec4(0.3, 0.1, 0.2, 0.3);

vec2 min2(vec2 a, vec2 b) {
    return a.x < b.x ? a : b;
}

float hash(const in vec3 p) {
    float h = dot(p, vec3(127.1, 311.7, 758.5453123));
    return fract(sin(h) * 43758.5453123);
}

// [iq] https://www.shadertoy.com/view/4sfGzS
float noiseText( in vec3 x) {
    vec3 p = floor(x), f = fract(x);
    f = f * f * (3. - f - f);
    vec2 uv = (p.xy + vec2(37., 17.) * p.z) + f.xy,
    //rg = textureLod( iChannel0, (uv+.5)/256., -100.).yx;
    rg = vec2(0.0, 0.0);
    return mix(rg.x, rg.y, f.z);
}

// ratio: ratio of hight/low frequencies
float fbmdust( in vec3 p, in float ratio) {
    return mix(noiseText(p * 3.), noiseText(p * 20.), ratio);
}

vec2 spiralArm( in vec3 p, in float thickness, in float blurAmout, in float blurStyle) {
    float dephase = 2.2, loop = 4.;
    float a = atan(p.x, p.z), // angle
    r = length(p.xz), lr = log(r), // distance to center
    th = (.1 - .25 * r), // thickness according to distance
    d = fract(.5 * (a - lr * loop) / pi); //apply rotation and scaling.
    d = (.5 / dephase - abs(d - .5)) * 2. * pi * r;
    d *= (1. - lr) / thickness; // space fct of distance
    // Perturb distance field
    float radialBlur = blurAmout * fbmdust(vec3(r * 4., 10. * d, 10. - 5. * p.y), blurStyle);
    return vec2(sqrt(d * d + 10. * p.y * p.y / thickness) - th * r * .2 - radialBlur);
}

vec2 dfGalaxy( in vec3 p, in float thickness, in float blurAmout, in float blurStyle) {
    return min2(spiralArm(p, thickness, blurAmout, blurStyle),
    spiralArm(vec3(p.z, p.y, -p.x), thickness, blurAmout, blurStyle));
}

vec2 map( in vec3 p) {
    float a = iMouse.x * .008 * pi + iTime * .3;
    p.xz = cos(a) * p.xz + sin(a) * vec2(p.xz.y, -p.xz.x);
    //    R(p.xz, iMouse.x*.008*pi+iTime*.3);
    return dfGalaxy(p, clamp(10. * sliderVal.x, .9, 10.), sliderVal.y, sliderVal.z);
}

//--------------------------------------------------------------

// assign color to the media
vec4 computeColor( in vec3 p, in float density, in float radius, in float id) {
    // color based on density alone, gives impression of occlusion within
    // the media
    vec4 result = mix(vec4(1., .9, .8, 1.), vec4(.4, .15, .1, 1.), density);
    // color added to the media
    result *= mix(colCenter,
    mix(colEdge2,
    mix(colEdge, colEdge3, step(.08, id)), step(-.05, id)),
    smoothstep(.2, .8, radius));
    return result;
}

// - Ray / Shapes Intersection -----------------------
bool sBox( in vec3 ro, in vec3 rd, in vec3 rad, out float tN, out float tF) {
    vec3 m = 1. / rd, n = m * ro,
    k = abs(m) * rad,
    t1 = -n - k, t2 = -n + k;
    tN = max(max(t1.x, t1.y), t1.z);
    tF = min(min(t2.x, t2.y), t2.z);
    return !(tN > tF || tF < 0.);
}

bool sSphere( in vec3 ro, in vec3 rd, in float r, out float tN, out float tF) {
    float b = dot(rd, ro), d = b * b - dot(ro, ro) + r;
    if (d < 0.) return false;
    tN = -b - sqrt(d);
    tF = -tN - b - b;
    return tF > 0.;
}

// ---------------------------------------------------
// Based on "Dusty nebula 4" by Duke (https://www.shadertoy.com/view/MsVXWW)
void mainImage(out vec4 fragColor, in vec2 fragCoord) {
    //    vec4 cSlider = processSliders(fragCoord);
    vec4 cSlider = vec4(0.0, 0.0, 0.0, 0.0);

    // camera
    float a = sliderVal.w * pi;
    vec3 ro = vec3(0., 2. * cos(a), -4.5 * sin(a)),
    ta = vec3(-.2, -.3, 0);

    // camera tx
    vec3 cw = normalize(ta - ro),
    cp = vec3(0., 1., 0.),
    cu = normalize(cross(cw, cp)),
    cv = normalize(cross(cu, cw));
    vec2 q = (fragCoord.xy) / iResolution.xy,
    p = -1. + 2. * q;
    p.x *= iResolution.x / iResolution.y;

    vec3 rd = normalize(p.x * cu + p.y * cv + 2.5 * cw);

    // ld, td: local, total density
    // w: weighting factor
    float ld = 0., td = 0., w = 0.;

    // t: length of the ray
    // d: distance function
    float d = 1., t = 0.;

    const float h = 0.1;

    vec4 sum = vec4(0);

    float min_dist = 0., max_dist = 0.,
    min_dist2 = 0., max_dist2 = 0.;

    if (sSphere(ro, rd, 4., min_dist, max_dist)) {
        if (sBox(ro, rd, vec3(4., 1.8, 4.), min_dist2, max_dist2)) {
            min_dist = max(.1, max(min_dist, min_dist2));
            max_dist = min(max_dist, max_dist2);

            t = min_dist * step(t, min_dist) + .1 * hash(rd + iTime);

            // raymarch loop
            vec4 col;
            for (int i = 0; i < 100; i++) {
                vec3 pos = ro + t * rd;

                // Loop break conditions.
                if (td > .9 || sum.a > .99 || t > max_dist) break;

                // evaluate distance function
                vec2 res = map(pos);
                d = max(res.x, .01);

                // point light calculations
                vec3 ldst = pos;
                ldst.y *= 1.6;
                vec3 ldst2 = pos;
                ldst2.y *= 3.6;
                float lDist = max(length(ldst), .1), //max(length(ldst), 0.001);
                lDist2 = max(length(ldst2), .1);
                // star in center
                vec3 lightColor = (1. - smoothstep(3., 4.5, lDist * lDist)) *
                mix(.015 * vec3(1., .5, .25) / (lDist * lDist),
                .02 * vec3(.5, .7, 1.) / (lDist2 * lDist2),
                smoothstep(.1, 2., lDist * lDist));
                sum.rgb += lightColor; //.015*lightColor/(lDist*lDist); // star itself and bloom around the light
                sum.a += .003 / (lDist * lDist);;

                if (d < h) {
                    // compute local density
                    ld = h - d;
                    // compute weighting factor
                    w = (1. - td) * ld;
                    // accumulate density
                    td += w + 1. / 60.;
                    // get color of object (with transparencies)
                    col = computeColor(pos, td, lDist * 2., res.y);
                    col.a *= td;
                    // colour by alpha
                    col.rgb *= col.a;
                    // alpha blend in contribution
                    sum += col * (1.0 - sum.a);
                }

                //float pitch = t/iResolution.x;
                //float dt = max(d * 0.25, .005); //pitch);
                // trying to optimize step size near the camera and near the light source
                t += max(d * .15 * max(min(length(ldst), length(ro)), 1.0), 0.005);
                td += .1 / 70.;
                //t += dt;
            }
            // simple scattering
            sum *= 1. / exp(ld * .2) * .8;
            sum = clamp(sum, 0., 1.);
        }
    }

    // Background color
    //    sum.rgb += vec3(clamp(2.*cos(.5*iTime),0.,.4))*(1. - sum.a)*pow(16.0*q.x*q.y*(1.-q.x)*(1.-q.y),.3);

    //Apply slider overlay

    float sum_f = length(sum) / 2.732; // sqrt(3)
    vec4 stars = clamp(starField(fragCoord), 0.0, 1.0) * (1.0 - sum_f);
    fragColor = vec4(sum.xyz + stars.xyz, 1.);
}

void main() {
    mainImage(gl_FragColor, gl_FragCoord.xy);
}