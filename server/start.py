from server.routes import *
from server.admin_routes import *

if __name__ == "__main__":
    import sklearn
    assert sklearn.__version__ == '1.1.0', ('Exact sklearn version is needed because PyTDC uses pickled random forest'
                                            'that requires exact version match.')

    with app.app_context():
        db.create_all()
    socketio.run(app, debug=True)
