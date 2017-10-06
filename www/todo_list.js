// create a binding object
var binding = new Shiny.InputBinding();

// add methods to it using jQuery's extend method
$.extend(binding, {

  find: function(scope) {

    // find all instances of class container
    return $(scope).find(".container");

  },

  // this method will be called on initialisation
  initialize: function(el){

     // extract the state from el
     // note here our bootstrapSwitch does not yet exist
     var state = $(el).data("state");

     // initialize our switch based on the extracted state
     // note $("#" + el.id) equals the input tag we generated
     $("#" + el.id).bootstrapSwitch("state",state);
     
     // now bootstrapSwitch does exist

  },

  // this method will also be called on initialisation (to pass the intial state to input$...)
  // and each time when the callback is triggered via the event bound in subscribe
  getValue: function(el) {

    // get the value from bootstrapSwitch
    var value = $(el).bootstrapSwitch('state');

    return value;
  },

  // we want to subscribe to the switchChange event
  // see http://bootstrapswitch.com/events.html
  subscribe: function(el, callback) {

    // only when the switchChange event is detected on instances of class bootstrapSwitch
    // trigger the getValue method and send the value to shiny
    $(document).on('switchChange.bootstrapSwitch', function(event){

      // callback which will tell Shiny to retrieve the value via getValue
      callback();
    });
  }
});

// register the binding so Shiny knows it exists
Shiny.inputBindings.register(binding);


var Task = function(name) {
    this.name = ko.observable(name);
}

var ViewModel = function() {
    var self = this;
    self.tasks = ko.observableArray([
        new Task("Get dog food"),
        new Task("Mow lawn"),
        new Task("Fix car"),
        new Task("Fix fence"),
        new Task("Walk dog"),
        new Task("Read book")
    ]);


    self.selectedTask = ko.observable();
    self.clearTask = function(data, event) {
        if (data === self.selectedTask()) {
            self.selectedTask(null);                
        }
        
        if (data.name() === "") {
           self.tasks.remove(data);   
        }
    };
    self.addTask = function() {
        var task = new Task("new");
        self.selectedTask(task);
        self.tasks.push(task);
    };

    self.isTaskSelected = function(task) {
       return task === self.selectedTask();  
    };
};

//control visibility, give element focus, and select the contents (in order)
ko.bindingHandlers.visibleAndSelect = {
    update: function(element, valueAccessor) {
        ko.bindingHandlers.visible.update(element, valueAccessor);
        if (valueAccessor()) {
            setTimeout(function() {
                $(element).find("input").focus().select();
            }, 0); //new tasks are not in DOM yet
        }
    }
};

ko.applyBindings(new ViewModel());