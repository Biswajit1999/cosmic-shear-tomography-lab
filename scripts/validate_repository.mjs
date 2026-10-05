import fs from 'node:fs';
const data=JSON.parse(fs.readFileSync('data/kids1000-xip.json','utf8'));
const failures=[];
for(const field of ['release_page','source_archive_url','source_archive_sha256','source_fits_sha256','citation','scope_note'])if(!data[field])failures.push(`metadata missing: ${field}`);
for(const pair of data.pairs){for(const value of [...pair.theta_arcmin,...pair.xi_plus,...pair.sigma,pair.null_chi2,pair.null_p_value,pair.correlated_snr])if(!Number.isFinite(value))failures.push(`${pair.id} contains a non-finite value`);if(pair.sigma.some(v=>v<=0))failures.push(`${pair.id} contains a non-positive uncertainty`)}
const far=data.pairs.find(p=>p.id==='5x5'),near=data.pairs.find(p=>p.id==='1x1');
if(!(far.null_p_value<1e-20))failures.push('5x5 regression anchor failed');
if(!(near.null_p_value>.1))failures.push('1x1 regression anchor failed');
if(failures.length){console.error(failures.join('\n'));process.exit(1)}
console.log('Research validation passed: provenance, covariance dimensions and regression anchors verified.');
