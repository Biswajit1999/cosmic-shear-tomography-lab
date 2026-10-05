const root=document.documentElement;
const button=document.getElementById('themeToggle');
const preferred=localStorage.getItem('cosmic-shear-theme')||(matchMedia('(prefers-color-scheme: dark)').matches?'dark':'light');
function apply(theme){root.dataset.theme=theme;button.textContent=theme==='dark'?'Day mode':'Night mode';button.setAttribute('aria-pressed',String(theme==='dark'));localStorage.setItem('cosmic-shear-theme',theme);window.dispatchEvent(new Event('lab:theme'))}
button.addEventListener('click',()=>apply(root.dataset.theme==='dark'?'light':'dark'));
apply(preferred);
