const express = require('express');
const bodyParser = require('body-parser');
const { exec } = require('child_process');
const path = require('path');

const app = express();
const port = 3000;


app.use(express.static(path.join(__dirname, 'public')));
app.use(bodyParser.json());


app.post('/calculate', (req, res) => {
    const { rho, mu, V, L } = req.body;
    const cmd = `matlab -batch "disp(interactiveCalculator(${rho}, ${mu}, ${V}, ${L}))"`

    exec(cmd, (error, stdout, stderr) => {
        if (error) {
            res.json({ error: stderr });
        } else {
            const lines = stdout.trim().split('\n');
            const result = parseFloat(lines[lines.length - 1]);
            res.json({ result });
            console.log("pass")
        }
    });
});

app.listen(port, () => {
    console.log(`Server running at http://localhost:${port}`);
});
