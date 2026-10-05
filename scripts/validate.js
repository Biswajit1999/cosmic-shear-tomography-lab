const fs=require('node:fs');
const crypto=require('node:crypto');
const required=['README.md','index.html','styles.css','app.js','data/kids1000-xip.json','research/source/kids_xipm.fits','scripts/derive_kids1000.py'];
const failures=[];
for(const file of required)if(!fs.existsSync(file))failures.push(`${file} missing`);
const data=JSON.parse(fs.readFileSync('data/kids1000-xip.json','utf8'));
if(data.pairs?.length!==15)failures.push('expected 15 tomographic pairs');
if(data.redshift_distributions?.length!==5)failures.push('expected five source redshift distributions');
for(const pair of data.pairs||[]){if(pair.theta_arcmin.length!==9||pair.xi_plus.length!==9||pair.sigma.length!==9)failures.push(`${pair.id} does not contain nine points`);if(pair.covariance.length!==9||pair.covariance.some(row=>row.length!==9))failures.push(`${pair.id} covariance is not 9x9`);if(!Number.isFinite(pair.null_p_value)||pair.null_p_value<0||pair.null_p_value>1)failures.push(`${pair.id} invalid p-value`)}
const hash=crypto.createHash('sha256').update(fs.readFileSync('research/source/kids_xipm.fits')).digest('hex');
if(hash!==data.source_fits_sha256)failures.push('FITS checksum does not match derivative metadata');
const combined=required.filter(f=>/\.(md|html|css|js|py)$/i.test(f)).map(f=>fs.readFileSync(f,'utf8')).join('\n').toLowerCase();
for(const token of ['to'+'do','place'+'holder','coming '+'soon'])if(combined.includes(token))failures.push(`unfinished token ${token}`);
if(failures.length){console.error(failures.join('\n'));process.exit(1)}
console.log(`KiDS-1000 validation passed: ${data.pairs.length} pairs, ${data.pairs.length*9} xi+ values, five n(z) distributions.`);
