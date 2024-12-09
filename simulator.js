const G = 6.6743 * 10 ** -11; // N kg^2 m^-2
const eta = 1; //softening factor
const density = 5000; // kg/m^3
const sphereScalar = (4 / 3) * Math.PI;
class Body {
    name = "";
    mass;
    radius;
    vX = 0;
    vY = 0;
    sX = 0;
    sY = 0;
    fX = 0;
    fY = 0;
    gradientX = 0;
    gradientY = 0;
    constructor(mass) {
        this.mass = mass;
        const volume = mass / density;
        this.radius = (volume / sphereScalar) ** (1 / 3);
    }
    timeStep(t) {
        const oldVX = (this.fX / this.mass) * t + this.vX; // acceleration
        const oldVY = (this.fY / this.mass) * t + this.vY;
        // const sX0 = this.sX;
        // const sY0 = this.sY;
        // this.sX += this.vX * t + (aX * t ** 2) / 2; // displacement
        // this.sY += this.vY * t + (aY * t ** 2) / 2;
        // this.vX += aX * t; // velocity
        // this.vY += aY * t;
        const Δx = this.getPosition(
            this.fX,
            this.gradientX,
            this.vX,
            this.sX,
            t
        );
        const Δy = this.getPosition(
            this.fY,
            this.gradientY,
            this.vY,
            this.sY,
            t
        );
        this.sX += Δx;
        this.sY += Δy;

        this.vX = this.getVelocity(
            this.fX,
            this.gradientX,
            this.mass,
            this.vX,
            Δx
        );

        this.vY = this.getVelocity(
            this.fY,
            this.gradientY,
            this.mass,
            this.vY,
            Δy
        );
    }
    getNetForce(bodies) {
        let fX = 0;
        let fY = 0;
        let gradX = 0;
        let gradY = 0;
        for (let b of bodies) {
            if (b === this) {
                continue;
            }
            const dX = b.sX - this.sX; //distance x
            const dY = b.sY - this.sY;
            const r = (dX ** 2 + dY ** 2) ** 0.5; //pythagorean theorem
            const f = this.gForce(this.mass, b.mass, r);
            const gradient = this.gForceGradient(this.mass, b.mass, r);
            const angle = Math.atan(Math.abs(dY) / Math.abs(dX));
            const sinAngle = Math.sin(angle);
            const cosAngle = Math.cos(angle);
            fX -= cosAngle * f * Math.sign(dX || 1);
            fY -= sinAngle * f * Math.sign(dY || 1);
            gradX += cosAngle * gradient * Math.sign(dX || 1);
            gradY += sinAngle * gradient * Math.sign(dY || 1);
        }
        this.fX = fX;
        this.fY = fY;
        this.gradientX = gradX;
        this.gradientY = gradY;

        return { fX, fY };
    }
    gForce(m1, m2, r) {
        return (-1 * (G * m1 * m2)) / (r ** 2 + eta ** 2);
    }
    gPotentialEnergy(m1, m2, r) {
        return -1 * G * m1 * m2 * (1 / r);
    }
    /**
     * @param {number} m1
     * @param {number} m2
     * @param {number} r
     * @return {number}
     */
    gForceGradient(m1, m2, r) {
        return (2 * G * m1 * m2) / r ** 3;
    }
    /**
     *
     *  @param {number} b b = f(x) - kx
     *  @param {number} k df/dx
     * @param {number} m mass
     * @param {number} u velocity
     * @param {number} Δs displacement
     * @returns {number}
     */
    getVelocity(b, k, m, u, Δs) {
        //work done
        const w = (k * Δs ** 2) / 2 + b * Δs;
        // ΔKE = 1/2 • m(v^2 - u^2)
        const v2 = u ** 2 + (w * 2) / m;
        const ΔsSign = Math.sign(Δs) || 1;
        const v2Sign = Math.sign(v2) || 1;
        if (this.name === "123") {
            console.log(
                this.name,
                "Δs",
                Δs,
                "b",
                b,
                "w",
                w,
                "v2",
                v2,
                "u",
                u,
                "v",
                Math.abs(v2) ** 0.5 * ΔsSign * v2Sign
            );
        }
        return Math.abs(v2) ** 0.5 * ΔsSign * v2Sign;
    }
    /**
     * @param {number} f - net force
     * @param {number} df_dx - force gradien - df/dx
     * @param {number} u - velocity
     * @param {number} x - displacement
     * @param {number} t - time
     * @returns {number}
     */
    getPosition(f, df_dx, u, x, t) {
        // F = ma
        // d2x/dt2 • m = kx + b
        // homogeneous solution:
        // xh = c1 • e^ωt + c2 • e^-ωt
        // where ω = √(k/m)
        // particular solution:
        // xp = -b/k
        // x = xh + xp
        // x0 = 0
        const k = df_dx;

        // approximate force as linear with respect to x
        // f = kx + b, where b is intercept and k is slope
        const b = f;
        // ω = iβ

        const β = Math.abs(k / this.mass) ** 0.5;
        // const c1 = 0.5 * (x + b / k + u / ω);
        // const c2 = 0.5 * (x + b / k - u / ω);

        // A = e^ωt, B = e^-ωt
        // 0.5((b/k + u/ω) e^ωt + (b/k - u/ω) e^-ωt)
        // 0.5((b/k) (e^ωt + e^-ωt) + (u/ω) (e^ωt - e^-ωt))
        // (b/k) cos(βt) + (u/ω) i sin(βt)
        // (b/k) cos(βt) + (u/β) sin(βt)

        // homogeneous solution to x
        const xh = (b / k) * Math.cos(β * t) + (u / β) * Math.sin(β * t);
        // particular solution
        const xp = (-1 * b) / k;
        return xh + xp;
    }
}

async function main() {
    // average mass of an astroid: 2E15 kg
    const b1 = new Body(2.1 * 10 ** 17);
    b1.name = "b1";
    const b2 = new Body(2.2 * 10 ** 14);
    b2.name = "b2";
    const b3 = new Body(2 * 10 ** 14);
    const b4 = new Body(2 * 10 ** 14);
    const b5 = new Body(2 * 10 ** 14);
    const b6 = new Body(2 * 10 ** 14);
    b2.sX = -500;
    b2.sY = 15;
    b2.vY = -100;
    b3.sY = -950;
    b3.sX = -550;
    b3.vY = -40;
    b4.sX = 1000;
    b4.vY = 100;
    b4.sY = -900;
    b5.sX = -700;
    b5.vY = 80;
    b5.sY = 1000;
    b6.sX = -900;
    b6.sY = -800;
    b6.vY = -70;

    const bodies = [b1, b2, b3, b4, b5, b6];
    //const bodies = [b1, b2, b3, b4];
    const chart = chartSetup();
    const config = getConfig();
    const chart2 = new Chart(document.getElementById("position"), config);
    chart2.options = {
        animation: {
            duration: 0
        }
    };

    displayBodies(chart, bodies);
    for (let i = 0; i < 10 ** 6; i++) {
        await timeStep(bodies);
        displayBodies(chart, bodies);
        speedDisplay(bodies);
        iterationDisplay(i);
        // let sum = 0;
        // bodies.forEach(b => {
        //     sum += b.pe + b.ke;
        // });
        // if (chart2.data.datasets[0].data.length > 500) {
        //     chart2.data.datasets[0].data.shift();
        // }
        // chart2.data.datasets[0].data.push({
        //     y: sum,
        //     x: i
        // });
        // chart2.update();
        // console.log(
        //     chart2.data.datasets[0].data[
        //         chart2.data.datasets[0].data.length - 1
        //     ]
        // );
    }
    console.log(bodies);
}
function removeOutOfBound(bodies, tolerance) {}
async function timeStep(bodies) {
    const time = 0.01;
    const waitTime = new Promise(resolve => {
        setTimeout(() => {
            resolve();
        }, time * 1000);
    });
    for (let b of bodies) {
        b.getNetForce(bodies);
    }
    for (let b of bodies) {
        b.timeStep(time);
    }

    return waitTime;
}
function iterationDisplay(i) {
    const div = document.getElementById("iteration-display");
    div.textContent = "iteration: " + i;
}
function speedDisplay(bodies) {
    const div = document.getElementById("speed-display");
    div.textContent = "";
    let i = 0;
    for (let b of bodies) {
        const p = document.createElement("p");
        p.textContent = `
        fX: ${b.fX.toExponential(1)},
        fY: ${b.fY.toExponential(1)},
        vX: ${b.vX.toExponential(1)},
        vY: ${b.vY.toExponential(1)},
        `;
        div.appendChild(p);
        i++;
    }
}
function displayBodies(chart, bodies) {
    try {
        // bodies.forEach((b, i) => {
        //     chart.options.plugins.annotation.annotations[i] = {
        //         radius: 0,
        //         type: "point",
        //         xMax: b.sX + b.radius,
        //         xMin: b.sX - b.radius,
        //         yMax: b.sY + b.radius,
        //         yMin: b.sY - b.radius,
        //         backgroundColor: "rgba(255, 99, 132, 0.25)"
        //     };
        // });
        chart.data.datasets = [
            {
                ...chart.data.datasets[0],
                data: bodies.map(b => ({
                    x: b.sX,
                    y: b.sY
                }))
            }
        ];

        chart.update();
    } catch (e) {
        console.error(e);
    }
}
function chartSetup() {
    const ctx = document.getElementById("myChart");
    return new Chart(ctx, getConfig());
}
function getConfig() {
    return {
        type: "scatter",
        data: getData(),
        options: {
            scales: {
                x: {
                    type: "linear",
                    max: 2000,
                    min: -2000,
                    position: "bottom"
                },
                y: {
                    type: "linear",
                    max: 2000,
                    min: -2000
                }
            },
            animation: {
                duration: 0
            },
            plugins: {
                annotation: {
                    annotations: {}
                }
            }
        }
    };
}
function getData() {
    return {
        datasets: [
            {
                label: "",
                data: [],
                backgroundColor: "rgb(255, 99, 132)"
            }
        ]
    };
}
try {
    main();
} catch (e) {
    console.error(e);
}
