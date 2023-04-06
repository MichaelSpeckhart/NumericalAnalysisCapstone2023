import React from 'react';
import { BrowserRouter, Route, Switch } from 'react-router-dom';
import {Home, Download, Demo, About} from './pages';

function App() {
  return (
    <BrowserRouter>
      <Switch>
        <Route exact path="/" component={Home} />
        <Route exact path="/download" component={Download} />
        <Route exact path="/demo" component={Demo} />
        <Route exact path="/about" component={About} />
      </Switch>
    </BrowserRouter>
  );
}

export default App;
