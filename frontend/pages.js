// npm install react-router-dom

import React from 'react';

function Home() {
  return (
    <div>
      <h1>Home Page</h1>
      <p>Welcome to the home page!</p>
    </div>
  );
}

function Download() {
    return (
      <div>
        <h1>Download Page</h1>
        <p>Welcome to the download page!</p>
      </div>
    );
}

function Demo() {
    fetch('example.wasm')
    .then(response => response.arrayBuffer())
    .then(buffer => WebAssembly.instantiate(buffer))
    .then(module => {
      // Use the WebAssembly module here
      const result = module.instance.exports.add(2, 3);
      console.log(result); // Output: 5
    });
}

function About() {
  return (
    <div>
      <h1>About Page</h1>
      <p>Learn more about us.</p>
    </div>
  );
}

export {Home, Download, Demo, About};
